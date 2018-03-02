#include <cassert>
#include <cmath>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "SphereSystem.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/IOHelper.hpp"

SphereSystem::SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv)
    : runConfig(configFile) {
    int mpiflag;
    MPI_Initialized(&mpiflag);
    assert(mpiflag);

    commRcp = getMPIWORLDTCOMM();

    stepCount = 0;
    snapID = 0;

    // TRNG pool must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(runConfig.rngSeed);
    collisionSolverPtr = std::make_shared<CollisionSolver>();
    collisionCollectorPtr = std::make_shared<CollisionCollector>();

    if (commRcp->getRank() == 0) {
        // prepare the output directory
        IOHelper::makeSubFolder("./result");
        // initialize on head rank
        setInitial(posFile);
    } else {
        sphere.reserve(runConfig.sphereNumber / commRcp->getSize() * 2);
    }

    // initial exchange of particles
    MPI_Barrier(MPI_COMM_WORLD);
    partition();
    MPI_Barrier(MPI_COMM_WORLD);

    printf("local: %d spheres on process %d\n", sphere.size(), commRcp->getRank());
    output();

    // initialize force calculate tree
    // initial resolve and reset move to zero
    // 5 pass to resolve initial overlaps

    printf("SphereSystem initialized on process: %d \n", commRcp->getRank());
}

void SphereSystem::setInitial(const std::string &initPosFile) {
    bool read = readXYZ(initPosFile);
    if (read) {
        return;
    }

    // this function executes only on process 0

    double minBoxEdge =
        std::min(runConfig.simBoxHigh[0] - runConfig.simBoxLow[0], runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]);
    minBoxEdge = std::min(runConfig.simBoxHigh[0] - runConfig.simBoxLow[0], minBoxEdge);

    if (runConfig.sphereRadiusSigmaHydro > 0) {
        rngPoolPtr->setLogNormalParameters(runConfig.sphereRadiusHydro, runConfig.sphereRadiusSigmaHydro);
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
    if (commRcp->getSize() > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &nSphere, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    if (commRcp->getRank() == 0) {
        fprintf(fptr, "%d\n", nSphere);
        fprintf(fptr, "SphereNumber: %d, time: %lf\n", nSphere, stepCount * runConfig.dt);
    }

    for (const auto &s : sphere) {
        fprintf(fptr, "S\t%8d\t%.6g\t%.6g\t%.6g\t%.6g\n", s.gid, s.radius, s.pos[0], s.pos[1], s.pos[2]);
    }

    fclose(fptr);
}

void SphereSystem::writeVTK(const std::string &baseFolder) {
    Sphere::writeVTP(sphere, baseFolder, std::to_string(snapID), commRcp->getRank());
    Sphere::writeVTU(sphere, baseFolder, std::to_string(snapID), commRcp->getRank());
    // write parallel head
    if (commRcp->getSize() > 1 && commRcp->getRank() == 0) {
        Sphere::writePVTP(baseFolder, "000", commRcp->getSize());
        std::vector<std::pair<int, std::string>> dataFields;
        std::vector<IOHelper::IOTYPE> types = {IOHelper::IOTYPE::Float64, IOHelper::IOTYPE::Float64};
        for (auto &layer : sphere[0].sphLayer) {
            dataFields.emplace_back(std::pair<int, std::string>(layer.second->dimension, layer.second->name));
        }
        dataFields.emplace_back(std::pair<int, std::string>(3, "stk"));
        Sphere::writePVTU(dataFields, types, baseFolder, std::to_string(snapID), commRcp->getSize());
    }
}

void SphereSystem::output() {
    // every 200 snaps in a sub folder
    const int num = 200;
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

Teuchos::RCP<TOP> SphereSystem::getMobOperator(bool manybody) const {
    Teuchos::RCP<TOP> mobOpRcp;
    if (manybody) {
        // get full mobility operator
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

Teuchos::RCP<TV> SphereSystem::getForceKnown() const {
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
            forcePtr(6 * i, c) = 0;
            forcePtr(6 * i + 1, c) = 0;
            forcePtr(6 * i + 2, c) = 0;
            // torque
            forcePtr(6 * i + 3, c) = 0;
            forcePtr(6 * i + 4, c) = 0;
            forcePtr(6 * i + 5, c) = 0;
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

        s.stepEuler(dt);
    }

    return;
}

void SphereSystem::resolveCollision(bool manybody, double buffer) {
    // positive buffer value means sphere collision radius is effectively smaller
    // i.e., less likely to collide

    // generate known velocity
    Teuchos::RCP<TV> forceKnownRcp = getForceKnown();
    Teuchos::RCP<TOP> mobOpRcp = getMobOperator(manybody && runConfig.hydro);
    Teuchos::RCP<TV> velocityKnownRcp = getVelocityKnown(mobOpRcp, forceKnownRcp);

    // Collect collision pair blocks
    std::vector<CollisionSphere> collisionSphereSrc;
    std::vector<CollisionSphere> collisionSphereTrg;
    auto &collector = *collisionCollectorPtr;
    collector.clear();

    interactManagerPtr->setupEssVec(collisionSphereSrc, collisionSphereTrg);
    printf("ESS vec created\n");
    auto nearInteractorPtr = interactManagerPtr->getNewNearInteraction();
    printf("nearInteractorPtr created\n");
    interactManagerPtr->setupNearInteractor(nearInteractorPtr, collisionSphereSrc, collisionSphereTrg);
    printf("setupNear\n");
    interactManagerPtr->calcNearInteraction(nearInteractorPtr, collisionSphereSrc, collisionSphereTrg, collector);
    printf("calcNear\n");

    // construct collision stepper
    collisionSolverPtr->setup(*(collector.collisionPoolPtr), sphereMobilityMapRcp, runConfig.dt, buffer);
    collisionSolverPtr->setControlLCP(1e-5, 200, false); // res, maxIte, NWTN refine
    collisionSolverPtr->solveCollision(mobOpRcp, velocityKnownRcp);

    return;
}

void SphereSystem::step() {
    resolveCollision(true, 0);
    // move forward
    Teuchos::RCP<TV> velocityRcp = Teuchos::rcp(new TV(*(collisionSolverPtr->getVelocityCol()), Teuchos::Copy));
    Teuchos::RCP<TV> velocityKnownRcp = collisionSolverPtr->getVelocityKnown();
    velocityRcp->update(1.0, *velocityKnownRcp, 1.0);
    moveEuler(velocityRcp);

    stepCount++;
    if (stepCount % runConfig.snapFreq == 0) {
        output();
    }
}