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

    // setup linear problem Ax=b
    AOpRcp = Teuchos::rcp(new SphereSTKSHOperator(spherePtr, name, fmmPtr, 0.5, 0, 1.0, 1.0));
    xRcp = Teuchos::rcp(new TMV(AOpRcp->getDomainMap(), 1, true));
    bRcp = Teuchos::rcp(new TMV(AOpRcp->getRangeMap(), 1, true));
    resRcp = Teuchos::rcp(new TMV(AOpRcp->getRangeMap(), 1, true));

    // testOperator();

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKMobMat::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const {
    // solve one linear problem for each column of X and Y
    const int nCol = X.getNumVectors();
    const int nRowLocal = X.getLocalLength();
    TEUCHOS_ASSERT(nCol == Y.getNumVectors());
    TEUCHOS_ASSERT(nRowLocal == Y.getLocalLength());

    // setup temporary space
    force.resize(mobMapRcp->getNodeNumElements());
    vel.resize(mobMapRcp->getNodeNumElements());

    TEUCHOS_ASSERT(nRowLocal == mobMapRcp->getNodeNumElements());
    TEUCHOS_ASSERT(nRowLocal == force.size());
    TEUCHOS_ASSERT(nRowLocal == vel.size());
    TEUCHOS_ASSERT(nRowLocal == 6 * spherePtr->size());

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();
    TEUCHOS_ASSERT(XPtr.dimension_0() == nRowLocal);
    TEUCHOS_ASSERT(YPtr.dimension_0() == nRowLocal);
    TEUCHOS_ASSERT(nCol == 1);

    for (int c = 0; c < nCol; c++) {
#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            force[i] = XPtr(i, c);
        }
        std::fill(vel.begin(), vel.end(), 0.0);

        solveMob(force.data(), vel.data());

#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            double temp = beta * YPtr(i, c) + alpha * vel[i];
            YPtr(i, c) = temp;
        }
    }
}

void SphereSTKMobMat::solveMob(const double *forcePtr, double *velPtr) const {

    if (commRcp->getRank() == 0)
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
    TEUCHOS_ASSERT(gridNorms.size() == 3 * gridWeights.size());
    TEUCHOS_ASSERT(gridCoords.size() == 3 * gridWeights.size());
    TEUCHOS_ASSERT(gridNumberIndex.size() == nLocal);
    TEUCHOS_ASSERT(gridNumberLength.size() == nLocal);
    TEUCHOS_ASSERT(sh.size() == nLocal);

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
            const int index = indexBase + j;
            // printf("indexBase %d, npts %d, index %d\n", indexBase, npts, index);
            Evec3 ro = A + B.cross(radius * Evec3(gridNorms[3 * (index)], gridNorms[3 * (index) + 1],
                                                  gridNorms[3 * (index) + 2]));
            rho[3 * (index) + 0] = ro[0];
            rho[3 * (index) + 1] = ro[1];
            rho[3 * (index) + 2] = ro[2];
        }
    }

    // step 2, compute - (1/2I+K) rho
    std::fill(b.begin(), b.end(), 0);
    AOpRcp->applyP2POP(rho.data(), b.data(), -0.5, 0, -1.0);
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
    // use b as temporary space
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBase = gridNumberIndex[i];
        const int npts = gridNumberLength[i];
        auto &density = sh[i].gridValue;
        TEUCHOS_ASSERT(npts * 3 == density.size());
        for (int j = 0; j < npts; j++) {
            const int index = indexBase + j;
            b[3 * (index) + 0] = density[3 * j + 0] - rho[3 * (index) + 0];
            b[3 * (index) + 1] = density[3 * j + 1] - rho[3 * (index) + 1];
            b[3 * (index) + 2] = density[3 * j + 2] - rho[3 * (index) + 2];
        }
    }

    AOpRcp->projectNullSpace(b.data());
    auto xPtr = xRcp->getLocalView<Kokkos::HostSpace>();
    xRcp->modify<Kokkos::HostSpace>();
    TEUCHOS_ASSERT(xPtr.dimension_0() == nRowLocal);
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        xPtr(i, 0) = b[i];
    }

    dumpTMV(bRcp, "B");
    dumpTMV(xRcp, "Xguess");
    // iterative solve
    // setup the problem
    Teuchos::RCP<Belos::LinearProblem<::TOP::scalar_type, TMV, TOP>> problemRcp =
        Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(AOpRcp, xRcp, bRcp));
    resRcp->putScalar(1.0);
    bool set = problemRcp->setProblem();
    TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error, "*** Belos::LinearProblem failed to set up correctly! ***");

    // testOperator();
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Num Blocks", 50); // larger values might trigger a std::bad_alloc inside Kokkos.
    solverParams->set("Maximum Iterations", 200);
    solverParams->set("Maximum Restarts", 50);
    solverParams->set("Convergence Tolerance", 1e-5);
    // solverParams->set("Orthogonalization", "ICGS");
    solverParams->set("Output Style", Belos::OutputType::General);
    // solverParams->set("Implicit Residual Scaling", "Norm of Initial Residual");
    // solverParams->set("Explicit Residual Scaling", "Norm of Initial Residual");
    // solverParams->set("Implicit Residual Scaling", "Norm of RHS");
    // solverParams->set("Explicit Residual Scaling", "Norm of RHS");
    solverParams->set("Implicit Residual Scaling", "None"); // default is preconditioned initial residual
    solverParams->set("Explicit Residual Scaling", "None");
    // solverParams->set("Output Frequency", 1);
    solverParams->set("Timer Label", "Iterative Mobility Solution");
    solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);

    Belos::SolverFactory<TOP::scalar_type, TMV, TOP> factory;
    Teuchos::RCP<Belos::SolverManager<TOP::scalar_type, TMV, TOP>> solverRcp = factory.create("GMRES", solverParams);
    solverRcp->setProblem(problemRcp);

    Belos::ReturnType result = solverRcp->solve();
    int numIters = solverRcp->getNumIters();
    if (commRcp->getRank() == 0) {
        std::cout << "Num of Iterations in Mobility Matrix: " << numIters << std::endl;
    }
    dumpTMV(xRcp, "Xsol");

    // step 3 compute velocity with solved density
    // rho+mu is saved in rho, the density generating the fluid velocity
#pragma omp parallel for
    for (int i = 0; i < nRowLocal; i++) {
        rho[i] += xPtr(i, 0);
        b[i] = 0;
    }

    // u is stored in b
    TEUCHOS_ASSERT(rho.size() == b.size());
    AOpRcp->applyP2POP(rho.data(), b.data(), 0, 1.0, 0);
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
            const int index = indexBase + j;
            Evec3 velpts(b[3 * index], b[3 * index + 1], b[3 * index + 2]);
            vel += gridWeights[index] * velpts;
            omega +=
                (radius * gridWeights[index]) *
                (Evec3(gridNorms[3 * (index)], gridNorms[3 * (index) + 1], gridNorms[3 * (index) + 2]).cross(velpts));
        }
        vel *= (1 / (4 * pi * radius * radius));
        omega *= (3 / (8 * pi * radius * radius * radius * radius));
        velPtr[6 * i + 0] = vel[0] * invVis;
        velPtr[6 * i + 1] = vel[1] * invVis;
        velPtr[6 * i + 2] = vel[2] * invVis;
        velPtr[6 * i + 3] = omega[0] * invVis;
        velPtr[6 * i + 4] = omega[1] * invVis;
        velPtr[6 * i + 5] = omega[2] * invVis;
    }
}

void SphereSTKMobMat::writeBackDensitySolution() {
    return;
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
            // rho must remain valid
            density[3 * j + 0] = rho[3 * (indexBase + j) + 0];
            density[3 * j + 1] = rho[3 * (indexBase + j) + 1];
            density[3 * j + 2] = rho[3 * (indexBase + j) + 2];
        }
    }
}

void SphereSTKMobMat::testOperator() {
    xRcp->putScalar(3.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x3b");
    xRcp->putScalar(1.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x1b");
    xRcp->putScalar(0.0);
    AOpRcp->apply(*xRcp, *bRcp);
    dumpTMV(bRcp, "x0b");
}