#include "SphereSTKMobMat.hpp"

constexpr double pi = 3.141592653589793238462643383279;

SphereSTKMobMat::SphereSTKMobMat(const std::vector<Sphere> *const spherePtr, const std::string name_,
                                 std::shared_ptr<STKFMM> &fmmPtr_, const double viscosity_)
    : spherePtr(spherePtr), fmmPtr(fmmPtr_), name(name_), viscosity(viscosity_) {
    // setup map
    commRcp = getMPIWORLDTCOMM();
    const int nLocal = spherePtr->size();

    sphereMapRcp = getTMAPFromLocalSize(nLocal, commRcp);
    mobMapRcp = getTMAPFromLocalSize(6 * nLocal, commRcp);

    // setup temporary space
    force.resize(mobMapRcp->getNodeNumElements(), 0);
    vel.resize(mobMapRcp->getNodeNumElements(), 0);

    // setup linear problem Ax=b
    AOpRcp = Teuchos::rcp(new SphereSTKSHOperator(spherePtr, name, fmmPtr, 0.5, 0, 0, 1.0, 1.0));
    xRcp = Teuchos::rcp(new TMV(AOpRcp->getDomainMap(), 1, true));
    bRcp = Teuchos::rcp(new TMV(AOpRcp->getRangeMap(), 1, true));

    // Create the GMRES solver.
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Num Blocks", 200); // larger than this might trigger a std::bad_alloc inside Kokkos.
    solverParams->set("Maximum Iterations", 1000);
    solverParams->set("Convergence Tolerance", 1e-6);
    solverParams->set("Timer Label", "Iterative Inverse Preconditioner");
    solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);

    solverRCP = factory.create("GMRES", solverParams);

    commRcp->barrier();

    // setup the problem
    problemRCP = Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(AOpRcp, xRcp, bRcp));
    problemRCP->setProblem();
    problemRCP->setRightPrec(precOp);
    solverRCP->setProblem(problemRCP);

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKMobMat::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const {
    // solve one linear problem for each column of X and Y
    const int nCol = X.getNumVectors();
    const int nRowLocal = X.getLocalLength();
    assert(nRowLocal == mobMapRcp->getNodeNumElements());

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();

    for (int c = 0; c < nCol; c++) {
#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            force[i] = XPtr(i, c);
        }

        solveMob(force.data(), vel.data());

#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            YPtr(i, c) = beta * YPtr(i, c) + alpha * vel[i];
        }
    }
}

void SphereSTKMobMat::solveMob(const double *forcePtr, double *velPtr) const {

    // step 1, compute rho
    const int nLocal = spherePtr->size();
    const auto &gridNorms = AOpRcp->getGridNorms();
    const auto &gridCoords = AOpRcp->getGridCoords();
    const auto &gridWeights = AOpRcp->getGridWeights();
    const auto &gridDofIndex = AOpRcp->getGridDofIndex();
    const auto &gridDofLength = AOpRcp->getGridDofLength();
    const auto &sh = AOpRcp->getSH();
    rho.clear();
    rho.resize(gridWeights.size() * 3);
    b.clear();
    b.resize(gridWeights.size() * 3);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        Evec3 f(forcePtr[6 * i], forcePtr[6 * i + 1], forcePtr[6 * i + 2]);     // force
        Evec3 t(forcePtr[6 * i + 3], forcePtr[6 * i + 4], forcePtr[6 * i + 5]); // torque
        const double radius = sh[i].radius;
        Evec3 A = f * (1 / (4 * pi * radius * radius));
        Evec3 B = t * (3 / (8 * pi * radius * radius * radius * radius));
        const int indexBase = gridDofIndex[i];
        const int npts = gridDofLength[i];
        for (int j = 0; j < npts; j++) {
            Evec3 ro = A + B.cross(radius * Evec3(gridNorms[3 * (indexBase + j)], gridNorms[3 * (indexBase + j) + 1],
                                                  gridNorms[3 * (indexBase + j) + 2]));
            rho[3 * (indexBase + j)] = ro[0];
            rho[3 * (indexBase + j) + 1] = ro[1];
            rho[3 * (indexBase + j) + 2] = ro[2];
        }
    }

    // step 2, compute - (1/2I+K) rho
    AOpRcp->runFMM(rho.data(), b.data(), -0.5, 0, 0, -1);
    auto bPtr = bRcp->getLocalView<Kokkos::HostSpace>();
    bRcp->modify<Kokkos::HostSpace>();
    const int nRowLocal = bRcp->getLocalLength();
    assert(nRowLocal == b.size());
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        bPtr(i, 0) = b[i];
    }

    // step 2 solve linear system
    solverRCP->reset(Belos::Problem);
    Belos::ReturnType result = solverRCP->solve();
    int numIters = solverRCP->getNumIters();
    if (commRcp->getRank() == 0) {
        std::cout << "Num of Iterations in Mobility Matrix: " << numIters << std::endl;
    }

    // step 3 compute velocity with solved density
    // rho+mu is saved in rho, the density generating the fluid velocity
    auto xPtr = xRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        rho[i] += xPtr(i, 0);
        b[i] = 0;
    }

    // u is stored in b
    AOpRcp->runFMM(rho.data(), b.data(), 0, 1, 0, 0);
    const double invVis = 1 / viscosity;
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBase = gridDofIndex[i];
        const int npts = gridDofLength[i];
        const double radius = sh[i].radius;
        Evec3 vel = Evec3::Zero();
        Evec3 omega = Evec3::Zero();
        for (int j = 0; j < npts; j++) {
            // velocity of this point
            Evec3 velpts(b[3 * (j + indexBase)], b[3 * (j + indexBase) + 1], b[3 * (j + indexBase) + 2]);
            vel += gridWeights[indexBase + j] * velpts;
            omega += radius * gridWeights[indexBase + j] *
                     (velpts.cross(Evec3(gridNorms[3 * (j + indexBase)], gridNorms[3 * (j + indexBase) + 1],
                                         gridNorms[3 * (j + indexBase) + 2])));
        }
        vel *= (1 / (4 * pi * radius * radius));
        omega *= -(1 / (4 * pi * radius * radius * radius * radius));
        velPtr[6 * indexBase] = vel[0] * invVis;
        velPtr[6 * indexBase + 1] = vel[1] * invVis;
        velPtr[6 * indexBase + 2] = vel[2] * invVis;
        velPtr[6 * indexBase + 3] = omega[0] * invVis;
        velPtr[6 * indexBase + 4] = omega[1] * invVis;
        velPtr[6 * indexBase + 5] = omega[2] * invVis;
    }
}