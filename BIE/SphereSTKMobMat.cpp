#include "SphereSTKMobMat.hpp"

constexpr double pi = 3.141592653589793238462643383279;

SphereSTKMobMat::SphereSTKMobMat(std::vector<Sphere> *const spherePtr, const std::string name_,
                                 std::shared_ptr<STKFMM> &fmmPtr_, const double viscosity_)
    : spherePtr(spherePtr), fmmPtr(fmmPtr_), name(name_), viscosity(viscosity_) {
    // setup map
    commRcp = getMPIWORLDTCOMM();
    const int nLocal = spherePtr->size();

    sphereMapRcp = getTMAPFromLocalSize(nLocal, commRcp);
    mobMapRcp = getTMAPFromLocalSize(6 * nLocal, commRcp);

    // setup temporary space
    force.resize(mobMapRcp->getNodeNumElements());
    vel.resize(mobMapRcp->getNodeNumElements());
    std::fill(force.begin(), force.end(), 0);
    std::fill(vel.begin(), vel.end(), 0);

    // setup linear problem Ax=b
    AOpRcp = Teuchos::rcp(new SphereSTKSHOperator(spherePtr, name, fmmPtr, 0.5, 0, 1.0, 1.0));
    xRcp = Teuchos::rcp(new TMV(AOpRcp->getDomainMap(), 1, true));
    bRcp = Teuchos::rcp(new TMV(AOpRcp->getRangeMap(), 1, true));

    testOperator();

    // setup the problem
    problemRcp = Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(AOpRcp, xRcp, bRcp));
    problemRcp->setProblem();

    // Create the GMRES solver.
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Num Blocks", 50); // larger than this might trigger a std::bad_alloc inside Kokkos.
    solverParams->set("Maximum Iterations", 1000);
    solverParams->set("Convergence Tolerance", 1e-6);
    solverParams->set("Timer Label", "Iterative Mobility Solution");
    solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);
    solverRcp = factory.create("GMRES", solverParams);
    solverRcp->setProblem(problemRcp);

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKMobMat::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const {
    // solve one linear problem for each column of X and Y
    const int nCol = X.getNumVectors();
    const int nRowLocal = X.getLocalLength();
    TEUCHOS_ASSERT(nRowLocal == mobMapRcp->getNodeNumElements());
    TEUCHOS_ASSERT(nRowLocal == force.size());
    TEUCHOS_ASSERT(nRowLocal == vel.size());
    TEUCHOS_ASSERT(nRowLocal == 6 * spherePtr->size());

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();

    for (int c = 0; c < nCol; c++) {
#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            force[i] = XPtr(i, c);
        }
        std::fill(vel.begin(), vel.end(), 0);

        solveMob(force.data(), vel.data());

#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            YPtr(i, c) = beta * YPtr(i, c) + alpha * vel[i];
        }
    }
}

void SphereSTKMobMat::solveMob(const double *forcePtr, double *velPtr) const {
    printf("start mobility solve\n");

    // step 1, compute rho
    const int nLocal = spherePtr->size();
    const auto &gridNorms = AOpRcp->getGridNorms();
    const auto &gridCoords = AOpRcp->getGridCoords();
    const auto &gridWeights = AOpRcp->getGridWeights();
    const auto &gridNumberIndex = AOpRcp->getGridNumberIndex();
    const auto &gridNumberLength = AOpRcp->getGridNumberLength();
    const auto &sh = AOpRcp->getSH();
    rho.resize(gridWeights.size() * 3);
    b.resize(gridWeights.size() * 3);
    std::fill(rho.begin(), rho.end(), 0);
    std::fill(b.begin(), b.end(), 0);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        Evec3 f(forcePtr[6 * i + 0], forcePtr[6 * i + 1], forcePtr[6 * i + 2]); // force
        Evec3 t(forcePtr[6 * i + 3], forcePtr[6 * i + 4], forcePtr[6 * i + 5]); // torque
        const double radius = sh[i].radius;
        Evec3 A = f * (1 / (4 * pi * radius * radius));
        Evec3 B = t * (3 / (8 * pi * radius * radius * radius * radius));
        const int indexBase = gridNumberIndex[i];
        const int npts = gridNumberLength[i];
        for (int j = 0; j < npts; j++) {
            Evec3 ro = A + B.cross(radius * Evec3(gridNorms[3 * (indexBase + j)], gridNorms[3 * (indexBase + j) + 1],
                                                  gridNorms[3 * (indexBase + j) + 2]));
            rho[3 * (indexBase + j) + 0] = ro[0];
            rho[3 * (indexBase + j) + 1] = ro[1];
            rho[3 * (indexBase + j) + 2] = ro[2];
        }
    }

    // step 2, compute - (1/2I+K) rho
    AOpRcp->applyP2POP(rho.data(), b.data(), -0.5, 0, -1);

    auto bPtr = bRcp->getLocalView<Kokkos::HostSpace>();
    bRcp->modify<Kokkos::HostSpace>();
    const int nRowLocal = bRcp->getLocalLength();
    TEUCHOS_ASSERT(nRowLocal == b.size());
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        bPtr(i, 0) = b[i];
    }

    // step 2 solve linear system
    // fill initial guess with current value in sh
    auto xPtr = xRcp->getLocalView<Kokkos::HostSpace>();
    xRcp->modify<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBase = gridNumberIndex[i];
        const int npts = gridNumberLength[i];
        auto &density = sh[i].gridValue;
        for (int j = 0; j < npts; j++) {
            const int index = indexBase + j;
            xPtr(3 * (index) + 0, 0) = density[3 * j + 0] - rho[3 * (index) + 0];
            xPtr(3 * (index) + 1, 0) = density[3 * j + 1] - rho[3 * (index) + 0];
            xPtr(3 * (index) + 2, 0) = density[3 * j + 2] - rho[3 * (index) + 0];
        }
    }
    dumpTMV(xRcp, "Xguess.mtx");

    // iterative solve

    solverRcp->reset(Belos::Problem);
    Belos::ReturnType result = solverRcp->solve();
    int numIters = solverRcp->getNumIters();
    if (commRcp->getRank() == 0) {
        std::cout << "Num of Iterations in Mobility Matrix: " << numIters << std::endl;
    }
    dumpTMV(bRcp, "B.mtx");
    dumpTMV(xRcp, "Xsol.mtx");

    // step 3 compute velocity with solved density
    // rho+mu is saved in rho, the density generating the fluid velocity
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        rho[i] += xPtr(i, 0);
        b[i] = 0;
    }

    // u is stored in b
    AOpRcp->applyP2POP(rho.data(), b.data(), 0, 1, 0);
    const double invVis = 1 / viscosity;
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBase = gridNumberIndex[i];
        const int npts = gridNumberLength[i];
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
        velPtr[6 * i + 0] = vel[0] * invVis;
        velPtr[6 * i + 1] = vel[1] * invVis;
        velPtr[6 * i + 2] = vel[2] * invVis;
        velPtr[6 * i + 3] = omega[0] * invVis;
        velPtr[6 * i + 4] = omega[1] * invVis;
        velPtr[6 * i + 5] = omega[2] * invVis;
    }
}

void SphereSTKMobMat::writeBackDensitySolution() {
    // write vector rho back to sphere
    auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
    const auto &gridNumberIndex = AOpRcp->getGridNumberIndex();
    const auto &gridNumberLength = AOpRcp->getGridNumberLength();

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBase = gridNumberIndex[i];
        const int npts = gridNumberLength[i];
        auto &density = sphere[i].getLayer(name).gridValue;
        TEUCHOS_ASSERT(density.size() == 3 * npts);
        for (int j = 0; j < npts; j++) {
            density[3 * j + 0] = rho[3 * (indexBase + j) + 0];
            density[3 * j + 1] = rho[3 * (indexBase + j) + 1];
            density[3 * j + 2] = rho[3 * (indexBase + j) + 2];
        }
    }
}

void SphereSTKMobMat::testOperator() {
    xRcp->putScalar(3.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x3b.mtx");
    xRcp->putScalar(1.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x1b.mtx");
    xRcp->putScalar(0.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x0b.mtx");
}