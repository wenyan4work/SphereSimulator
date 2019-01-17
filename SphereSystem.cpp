#include "SphereSystem.hpp"

#include <cassert>
#include <cmath>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "BIE/SphereSTKMobMat.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/IOHelper.hpp"

SphereSystem::SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv)
    : runConfig(configFile) {
    int mpiflag;
    MPI_Initialized(&mpiflag);
    assert(mpiflag);

    Eigen::initParallel();
    Eigen::setNbThreads(1);   // disable threading in eigen
    srand(runConfig.rngSeed); // Eigen uses srand internally

    commRcp = getMPIWORLDTCOMM();

    stepCount = 0;
    snapID = 0;

    // TRNG pool & fmm must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(runConfig.rngSeed);
    collisionSolverPtr = std::make_shared<CollisionSolver>();
    collisionCollectorPtr = std::make_shared<CollisionCollector>();

    if (commRcp->getRank() == 0) {
        // prepare the output directory
        IOHelper::makeSubFolder("./result");
        // initialize on head rank
        setInitial(posFile);
        // initialize output data fields
        // used only on rank 0

    } else {
        sphere.reserve(runConfig.sphereNumber / commRcp->getSize() * 2);
    }

    // initial exchange of particles
    commRcp->barrier();
    partition();
    commRcp->barrier();

    // activate PVel and Traction operator
    if (runConfig.hydro) {
        fmmPtr = std::make_shared<stkfmm::STKFMM>(runConfig.pFMM, 4000, stkfmm::PAXIS::NONE, 9);
        for (auto &s : sphere) {
            s.addLayer("stkmob", Shexp::KIND::STK, runConfig.pSH, s.radius, Equatn::UnitRandom());
            s.addLayer("stkcol", Shexp::KIND::STK, runConfig.pSH, s.radius, Equatn::UnitRandom());
        }
        // manually initialize it on all ranks
        // VTU data always in Float64 format
        dataFieldVTU.emplace_back(3, IOHelper::IOTYPE::Float64, "stkmob");
        dataFieldVTU.emplace_back(3, IOHelper::IOTYPE::Float64, "stkcol");
    }

    // initialize force calculate tree
    // initial resolve and reset move to zero
    applyMonoLayer();

    printf("local: %lu spheres on process %d\n", sphere.size(), commRcp->getRank());
    output();

    printf("SphereSystem initialized on process: %d \n", commRcp->getRank());
}

void SphereSystem::applyMonoLayer() {
    if (!runConfig.monolayer) {
        return;
    }
    const int nLocal = sphere.size();
    const double monoZ = 0.5 * (runConfig.simBoxLow[2] + runConfig.simBoxHigh[2]);
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sphere[i].pos[2] = monoZ;
    }
}

void SphereSystem::setInitial(const std::string &initPosFile) {
    bool read = readXYZ(initPosFile);
    if (read) {
        return;
    }

    // this function executes only on process 0

    double minBoxEdge =
        std::min(runConfig.simBoxHigh[0] - runConfig.simBoxLow[0], runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]);
    minBoxEdge = std::min(runConfig.simBoxHigh[2] - runConfig.simBoxLow[2], minBoxEdge);

    if (runConfig.sphereRadiusSigmaHydro > 0) {
        rngPoolPtr->setLogNormalParameters(log(runConfig.sphereRadiusHydro), runConfig.sphereRadiusSigmaHydro);
    }

    const int nSphereGlobal = runConfig.sphereNumber;

    for (int i = 0; i < nSphereGlobal; i++) {
        double radius = runConfig.sphereRadiusSigmaHydro > 0 ? rngPoolPtr->getLN(0) : runConfig.sphereRadiusHydro;
        double px = rngPoolPtr->getU01(0) * (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0]) + runConfig.simBoxLow[0];
        double py = rngPoolPtr->getU01(0) * (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) + runConfig.simBoxLow[1];
        double pz = rngPoolPtr->getU01(0) * (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]) + runConfig.simBoxLow[2];
        Equatn orientation;
        EquatnHelper::setUnitRandomEquatn(orientation, rngPoolPtr->getU01(0), rngPoolPtr->getU01(0),
                                          rngPoolPtr->getU01(0));
        sphere.emplace_back(i, radius, radius * runConfig.sphereRadiusCollisionRatio, Evec3(px, py, pz), orientation);
    }

    // check volume fraction
    double particleVolume = 0;
    for (int i = 0; i < nSphereGlobal; i++) {
        const auto &s = sphere[i];
        particleVolume += 3.1416 * (4 / 3.0) * pow(s.radius, 3);
    }
    double boxVolume = (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0]) *
                       (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) *
                       (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]);

    std::cout << "initial volume fraction: " << particleVolume / boxVolume << std::endl;
}

bool SphereSystem::readXYZ(const std::string &filename) {
    std::ifstream myfile(filename);
    if (!myfile.is_open()) {
        printf("Error: cannot open initial position file. use random initial configuration\n");
        return false;
    }
    std::string line;
    getline(myfile, line); //
    getline(myfile, line); // read two header lines

    while (std::getline(myfile, line)) {
        // read 5 numbers: gid, radius, x,y,z
        // set radiusCollision = radius
        std::istringstream liness(line);
        int gid;
        char type;
        double px, py, pz;
        double radius;
        liness >> type >> gid >> radius >> px >> py >> pz;
        std::cout << gid << " " << radius << " " << px << " " << py << " "
                  << " " << pz << "\n";
        sphere.emplace_back(gid, radius, radius * runConfig.sphereRadiusCollisionRatio, Evec3(px, py, pz));
    }
    myfile.close();

    std::cout << "read initial finished" << std::endl;
    return true;
}

void SphereSystem::writeXYZ(const std::string &baseFolder) {
    // each rank write a xyz file separately. rank 0 write the two lines of header
    std::string name = baseFolder + std::string("SphereXYZ_") + std::string("r") + std::to_string(commRcp->getRank()) +
                       "_" + std::to_string(snapID) + ".xyz";
    FILE *fptr = fopen(name.c_str(), "w");
    const int nSphereLocal = sphere.size();
    int nSphere = nSphereLocal;
    // if (commRcp->getSize() > 1) {
    // MPI_Allreduce(MPI_IN_PLACE, &nSphere, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // }

    if (commRcp->getRank() == 0) {
        fprintf(fptr, "%d\n", sphereMapRcp->getGlobalNumElements());
        fprintf(fptr, "SphereNumber: %d, time: %lf\n", nSphere, stepCount * runConfig.dt);
    }

    for (const auto &s : sphere) {
        fprintf(fptr, "S\t%d\t%.14g\t%.14g\t%.14g\t%.14g\n", s.gid, s.radius, s.pos[0], s.pos[1], s.pos[2]);
    }

    fclose(fptr);
}

void SphereSystem::writeVTK(const std::string &baseFolder) {
    // write parallel head
    if (commRcp->getRank() == 0) {
        Sphere::writePVTP(baseFolder, std::to_string(snapID), commRcp->getSize());
        // dataFields on all ranks should have been setup during initialization
        Sphere::writePVTU(dataFieldVTU, baseFolder, std::to_string(snapID), commRcp->getSize());
        CollisionCollector::writePVTP(baseFolder, std::to_string(snapID), commRcp->getSize());
    }

    // write files from each rank
    Sphere::writeVTP(sphere, baseFolder, std::to_string(snapID), commRcp->getRank());
    Sphere::writeVTU(sphere, dataFieldVTU, baseFolder, std::to_string(snapID), commRcp->getRank());
    CollisionCollector::writeVTP(*(collisionCollectorPtr->collisionPoolPtr), baseFolder, std::to_string(snapID),
                                 commRcp->getRank());
}

void SphereSystem::output() {
    const int num = std::max(200 / commRcp->getSize(), 1); // limit max number of files per folder
    int k = snapID / num;
    int low = k * num, high = k * num + num - 1;
    std::string baseFolder =
        "./result/result" + std::to_string(low) + std::string("-") + std::to_string(high) + std::string("/");
    IOHelper::makeSubFolder(baseFolder);
    writeXYZ(baseFolder);
    writeVTK(baseFolder);
    snapID++;
}

void SphereSystem::partition() {
    // initialize
    if (!interactManagerPtr) {
        interactManagerPtr = std::make_shared<InteractionManager<double, 3, Sphere, Sphere>>(&sphere, &sphere);
    }

    // create a nearInteractor for full particle
    auto nearInteractFullParPtr = interactManagerPtr->getNewNearInteraction();
    interactManagerPtr->partitionObject(nearInteractFullParPtr);

    // setup the new sphereMap
    sphereMapRcp = getTMAPFromLocalSize(sphere.size(), commRcp);
    sphereMobilityMapRcp = getTMAPFromLocalSize(sphere.size() * 6, commRcp);

    // setup the globalIndex
    int globalIndexBase = sphereMapRcp->getMinGlobalIndex(); // this is a contiguous map
    const int nsphere = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < nsphere; i++) {
        sphere[i].globalIndex = i + globalIndexBase;
    }

    return;
}

Teuchos::RCP<TOP> SphereSystem::getMobOperator(bool manybody, std::string name) {
    Teuchos::RCP<TOP> mobOpRcp;
    if (manybody) {
        // get full mobility operator
        mobOpRcp = Teuchos::rcp(new SphereSTKMobMat(&sphere, name, fmmPtr, runConfig.viscosity));
    } else {
        Teuchos::RCP<TCMAT> mobMatRcp;
        // get mobility matrix for local drag only
        const double mu = runConfig.viscosity;
        const int nSphereLocal = sphere.size();

        // diagonal hydro mobility operator, with rotation
        const int localSize = nSphereLocal * 6;
        Kokkos::View<size_t *> rowPointers("rowPointers", localSize + 1);
        rowPointers[0] = 0;
        for (int i = 1; i <= localSize; i++) {
            rowPointers[i] = rowPointers[i - 1] + 1;
        }
        Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localSize]);
        Kokkos::View<double *> values("values", rowPointers[localSize]);
        for (int i = 0; i < rowPointers[localSize]; i++) {
            columnIndices[i] = i;
        }

#pragma omp parallel for
        for (int i = 0; i < nSphereLocal; i++) {
            const auto &s = sphere[i];
            const double radius = s.radius;
            const double invDragTrans = 1.0 / (6 * Pi * radius * mu);
            const double invDragRot = 1.0 / (8 * Pi * radius * radius * radius * mu);
            values[6 * i] = invDragTrans;
            values[6 * i + 1] = invDragTrans;
            values[6 * i + 2] = invDragTrans;
            values[6 * i + 3] = invDragRot;
            values[6 * i + 4] = invDragRot;
            values[6 * i + 5] = invDragRot;
        }

        mobMatRcp =
            Teuchos::rcp(new TCMAT(sphereMobilityMapRcp, sphereMobilityMapRcp, rowPointers, columnIndices, values));
        mobMatRcp->fillComplete(sphereMobilityMapRcp, sphereMobilityMapRcp); // domainMap, rangeMa
        mobOpRcp = mobMatRcp;
    }

    return mobOpRcp;
}

Teuchos::RCP<TV> SphereSystem::getVelocityKnown(Teuchos::RCP<TOP> &mobilityOpRcp, Teuchos::RCP<TV> &forceRcp) const {
    // velocity from known forcing
    Teuchos::RCP<TV> velocityKnownRcp = Teuchos::rcp<TV>(new TV(sphereMobilityMapRcp, true));
    mobilityOpRcp->apply(*forceRcp, *velocityKnownRcp);

    // extra velocity, e.g., Brownian noise
    if (runConfig.scaleBrown > 0) {
        Teuchos::RCP<TV> velocityBrownRcp = getVelocityBrown();
        velocityKnownRcp->update(1.0, *velocityBrownRcp, 1.0);
    }

    return velocityKnownRcp;
}

Teuchos::RCP<TV> SphereSystem::getForceKnown() {
    // 6 dof per sphere, force+torque
    Teuchos::RCP<TV> forceKnownRcp = Teuchos::rcp<TV>(new TV(sphereMobilityMapRcp, true));

    auto forcePtr = forceKnownRcp->getLocalView<Kokkos::HostSpace>();
    forceKnownRcp->modify<Kokkos::HostSpace>();

    const int sphereLocalNumber = sphere.size();
    assert(forcePtr.dimension_0() == sphereLocalNumber * 6);
    assert(forcePtr.dimension_1() == 1);

    for (int c = 0; c < forcePtr.dimension_1(); c++) {
#pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < sphereLocalNumber; i++) {
            // force
            sphere[i].forceNonCol[0] = (forcePtr(6 * i, c) = runConfig.extForce[0]);
            sphere[i].forceNonCol[1] = (forcePtr(6 * i + 1, c) = runConfig.extForce[1]);
            sphere[i].forceNonCol[2] = (forcePtr(6 * i + 2, c) = runConfig.extForce[2]);
            // torque
            sphere[i].torqueNonCol[0] = (forcePtr(6 * i + 3, c) = runConfig.extTorque[0]);
            sphere[i].torqueNonCol[1] = (forcePtr(6 * i + 4, c) = runConfig.extTorque[1]);
            sphere[i].torqueNonCol[2] = (forcePtr(6 * i + 5, c) = runConfig.extTorque[2]);
        }
    }

    commRcp->barrier();

    return forceKnownRcp;
}

Teuchos::RCP<TV> SphereSystem::getVelocityBrown() const {
    Teuchos::RCP<TV> velocityBrownRcp = Teuchos::rcp<TV>(new TV(sphereMobilityMapRcp, false));

    // extra velocity, e.g., Brownian noise
    auto velocityPtr = velocityBrownRcp->getLocalView<Kokkos::HostSpace>();
    velocityBrownRcp->modify<Kokkos::HostSpace>();

    const int sphereLocalNumber = sphere.size();
    assert(velocityPtr.dimension_0() == sphereLocalNumber * 6);
    assert(velocityPtr.dimension_1() == 1);

    const double mu = runConfig.viscosity;
    const double fackBT = sqrt(2 * runConfig.kBT / runConfig.dt) * runConfig.scaleBrown;

    for (int c = 0; c < velocityPtr.dimension_1(); c++) {
#pragma omp parallel for schedule(dynamic, 1024)
        for (int i = 0; i < sphereLocalNumber; i++) {
            const int threadId = omp_get_thread_num();
            const auto &s = sphere[i];
            const double radius = s.radius;
            const double invDragTrans = 1.0 / (6 * Pi * radius * mu);
            const double invDragRot = 1.0 / (8 * Pi * radius * radius * radius * mu);
            const double dx = sqrt(invDragTrans) * fackBT;
            const double dr = sqrt(invDragRot) * fackBT;
            // translation
            velocityPtr(6 * i, c) = dx * rngPoolPtr->getN01(threadId);
            velocityPtr(6 * i + 1, c) = dx * rngPoolPtr->getN01(threadId);
            velocityPtr(6 * i + 2, c) = dx * rngPoolPtr->getN01(threadId);
            // rotation
            velocityPtr(6 * i + 3, c) = dr * rngPoolPtr->getN01(threadId);
            velocityPtr(6 * i + 4, c) = dr * rngPoolPtr->getN01(threadId);
            velocityPtr(6 * i + 5, c) = dr * rngPoolPtr->getN01(threadId);
        }
    }

    return velocityBrownRcp;
}

void SphereSystem::moveEuler(Teuchos::RCP<TV> &velocityRcp) {

    assert(velocityRcp->getMap()->getNodeNumElements() == sphere.size() * 6);
    auto velocityPtr = velocityRcp->getLocalView<Kokkos::HostSpace>();
    velocityRcp->modify<Kokkos::HostSpace>();

    const int sphereLocalNumber = sphere.size();
    assert(velocityPtr.dimension_0() == sphereLocalNumber * 6);
    assert(velocityPtr.dimension_1() == 1);

    const int c = 0; // only 1 column in the TV
    const double dt = runConfig.dt;

#pragma omp parallel for schedule(dynamic, 1024)
    for (int i = 0; i < sphereLocalNumber; i++) {
        // translation
        const auto vx = velocityPtr(6 * i, c);
        const auto vy = velocityPtr(6 * i + 1, c);
        const auto vz = velocityPtr(6 * i + 2, c);
        // rotation
        const auto wx = velocityPtr(6 * i + 3, c);
        const auto wy = velocityPtr(6 * i + 4, c);
        const auto wz = velocityPtr(6 * i + 5, c);
        // update
        auto &s = sphere[i];
        s.vel = Evec3(vx, vy, vz);
        s.omega = Evec3(wx, wy, wz);
        // printf("%d vel: %.14g,%.14g,%.14g, omega: %.14g, %.14g, %.14g\n", i, vx, vy, vz, wx, wy, wz);
        s.stepEuler(dt);
    }

    applyMonoLayer();
}

void SphereSystem::resolveCollision(bool manybody, double buffer) {
    // positive buffer value means sphere collision radius is effectively smaller
    // i.e., less likely to collide
    if (fmmPtr) {
        fitFMMBox();
    }

    // generate known velocity
    Teuchos::RCP<TV> forceKnownRcp = getForceKnown();
    Teuchos::RCP<TOP> mobOpRcp = getMobOperator(manybody && runConfig.hydro, std::string("stkmob"));
    Teuchos::RCP<TV> velocityKnownRcp = getVelocityKnown(mobOpRcp, forceKnownRcp);
    if (runConfig.hydro && manybody) {
        Teuchos::RCP<SphereSTKMobMat> stkmobopRcp = rcp_dynamic_cast<SphereSTKMobMat>(mobOpRcp, true);
        stkmobopRcp->writeBackDensitySolution();
    }

    // Collect collision pair blocks
    auto &collector = *collisionCollectorPtr;
    collector.clear();

    // no collision for only 1 object
    // temporary workaround for a bug in interact manager
    std::vector<CollisionSphere> collisionSphereSrc;
    std::vector<CollisionSphere> collisionSphereTrg;

    interactManagerPtr->setupEssVec(collisionSphereSrc, collisionSphereTrg);
    if (commRcp->getRank() == 0)
        printf("ESS vec created\n");

    auto nearInteractorPtr = interactManagerPtr->getNewNearInteraction();
    if (commRcp->getRank() == 0)
        printf("nearInteractorPtr created\n");

    interactManagerPtr->setupNearInteractor(nearInteractorPtr, collisionSphereSrc, collisionSphereTrg);
    if (commRcp->getRank() == 0)
        printf("setupNear\n");

    interactManagerPtr->calcNearInteraction(nearInteractorPtr, collisionSphereSrc, collisionSphereTrg, collector);
    if (commRcp->getRank() == 0)
        printf("calcNear\n");

    // construct collision stepper
    collisionSolverPtr->setup(*(collector.collisionPoolPtr), sphereMobilityMapRcp, runConfig.dt, buffer);
    collisionSolverPtr->setControlLCP(1e-5, 2000, false); // res, maxIte, NWTN refine

    mobOpRcp = getMobOperator(manybody && runConfig.hydro, std::string("stkcol"));
    collisionSolverPtr->solveCollision(mobOpRcp, velocityKnownRcp);
    collisionSolverPtr->writebackGamma(*(collisionCollectorPtr->collisionPoolPtr));
    if (runConfig.hydro && manybody) {
        Teuchos::RCP<SphereSTKMobMat> stkmobopRcp = rcp_dynamic_cast<SphereSTKMobMat>(mobOpRcp, true);
        stkmobopRcp->writeBackDensitySolution();
    }

    // writeBack force torque solution
    const int nLocal = sphere.size();
    auto forceColRcp = collisionSolverPtr->getForceCol();

    auto forceColPtr = forceColRcp->getLocalView<Kokkos::HostSpace>();
    forceColRcp->modify<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sphere[i].forceCol[0] = forceColPtr(6 * i + 0, 0);
        sphere[i].forceCol[1] = forceColPtr(6 * i + 1, 0);
        sphere[i].forceCol[2] = forceColPtr(6 * i + 2, 0);
        sphere[i].torqueCol[0] = forceColPtr(6 * i + 3, 0);
        sphere[i].torqueCol[1] = forceColPtr(6 * i + 4, 0);
        sphere[i].torqueCol[2] = forceColPtr(6 * i + 5, 0);
    }

    return;
}

void SphereSystem::step() {
    if (commRcp->getRank() == 0)
        std::cout << "RECORD: TIMESTEP " << stepCount << " TIME " << stepCount * runConfig.dt << std::endl;

    resolveCollision(true, 0);

    stepCount++;
    if (stepCount % runConfig.snapFreq == 0) {
        output();
    }

    // move forward
    Teuchos::RCP<TV> velocityRcp = Teuchos::rcp(new TV(*(collisionSolverPtr->getVelocityCol()), Teuchos::Copy));
    Teuchos::RCP<TV> velocityKnownRcp = collisionSolverPtr->getVelocityKnown();
    velocityRcp->update(1.0, *velocityKnownRcp, 1.0);
    moveEuler(velocityRcp);
}

void SphereSystem::calcBoundaryCollision() {
    // a demo of how to calculate boundary collisions
    auto collisionPoolPtr = collisionCollectorPtr->collisionPoolPtr; // shared_ptr
    const int nThreads = collisionPoolPtr->size();
    const int nLocal = sphere.size();

    const int maxGlobalIndexOnLocal = sphereMapRcp->getMaxGlobalIndex();
    const int minGlobalIndexOnLocal = sphereMapRcp->getMinGlobalIndex();

    // a spherical shell boundary
    const Evec3 shellCenter(runConfig.simBoxHigh[0] * 0.5, runConfig.simBoxHigh[1] * 0.5,
                            runConfig.simBoxHigh[2] * 0.5);
    const double shellRadius = runConfig.simBoxHigh[2];

#pragma omp parallel num_threads(nThreads)
    {
        const int threadId = omp_get_thread_num();

#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &s = sphere[i];

            // do this for each boundary. add as many boundaries as you want
            {
                // calculate collision location
                Evec3 rvec = s.pos - shellCenter;
                double rnorm = rvec.norm();
                if (rnorm > shellRadius - s.radiusCollision) {
                    Evec3 normI = (-rvec).normalized();
                    double phi0 = -(rnorm - shellRadius); // negative
                    double gammaGuess = -phi0;            // positive
                    // add a new collision block. this block has only 6 non zero entries.
                    // passing sy.gid+1/globalIndex+1 as a 'fake' colliding body j, which is actually not used in the
                    // solver when oneside=true, out of range index is ignored
                    (*collisionPoolPtr)[threadId].emplace_back(phi0, -phi0, s.gid, s.gid + 1, s.globalIndex,
                                                               s.globalIndex + 1, normI, Evec3(0, 0, 0),
                                                               normI * s.radiusCollision, Evec3(0, 0, 0), true);
                }
            }
        }
    }
}

void SphereSystem::fitFMMBox() {
    if (runConfig.xPeriodic * runConfig.yPeriodic * runConfig.zPeriodic != 0) {
        // if any direction is periodic, use the value set in runConfig.txt
        fmmPtr->setBox(runConfig.simBoxLow[0], runConfig.simBoxHigh[0], runConfig.simBoxLow[1], runConfig.simBoxHigh[1],
                       runConfig.simBoxLow[2], runConfig.simBoxHigh[2]);
        return;
    }

    // otherwise, automatically fit the periodic box
    double xlow, ylow, zlow;
    double xhigh, yhigh, zhigh;
    xlow = ylow = zlow = 1e9;
    xhigh = yhigh = zhigh = -1e9;

    const int nLocal = sphere.size();
#pragma omp parallel for reduction(min : xlow, ylow, zlow) reduction(max : xhigh, yhigh, zhigh)
    for (int i = 0; i < nLocal; i++) {
        const auto &s = sphere[i];
        xlow = std::min(xlow, s.pos[0] - s.radius * 2);
        ylow = std::min(ylow, s.pos[1] - s.radius * 2);
        zlow = std::min(zlow, s.pos[2] - s.radius * 2);
        xhigh = std::max(xhigh, s.pos[0] + s.radius * 2);
        yhigh = std::max(yhigh, s.pos[1] + s.radius * 2);
        zhigh = std::max(zhigh, s.pos[2] + s.radius * 2);
    }

    // global min/max
    double localLow[3] = {xlow, ylow, zlow};
    double localHigh[3] = {xhigh, yhigh, zhigh};
    double globalLow[3] = {xlow, ylow, zlow};
    double globalHigh[3] = {xhigh, yhigh, zhigh};

    // add some buffer around the box
    for (int k = 0; k < 3; k++) {
        const double buffer = fabs(0.1 * (localHigh[k] - localLow[k]));
        localLow[k] -= buffer;
        localHigh[k] += buffer;
    }

    // printf("Local box bound low: %lf,%lf,%lf\n", localLow[0], localLow[1], localLow[2]);
    // printf("Local box bound high: %lf,%lf,%lf\n", localHigh[0], localHigh[1], localHigh[2]);

    Teuchos::reduceAll(*commRcp, Teuchos::MinValueReductionOp<int, double>(), 3, localLow, globalLow);
    Teuchos::reduceAll(*commRcp, Teuchos::MaxValueReductionOp<int, double>(), 3, localHigh, globalHigh);

    if (commRcp->getRank() == 0) {
        printf("Global box bound low: %lf,%lf,%lf\n", globalLow[0], globalLow[1], globalLow[2]);
        printf("Global box bound high: %lf,%lf,%lf\n", globalHigh[0], globalHigh[1], globalHigh[2]);
    }

    fmmPtr->setBox(globalLow[0], globalHigh[0], globalLow[1], globalHigh[1], globalLow[2], globalHigh[2]);
}