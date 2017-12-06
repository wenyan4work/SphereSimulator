//#pragma GCC push_options
//// for debug
//#pragma GCC optimize ("O0")

#include "SphereSystem.hpp"

#include "CPSolver.hpp"
#include "GeoUtil.hpp"
#include "TpetraUtil.hpp"
#include "ZDD.hpp"

#include <cmath>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <sys/stat.h>
#include <unordered_map>

#include <random>

#include "mpi.h"

template<typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

using namespace HydroSphere;

void periodicWrap(double &pos, const double &trg, const double &pbcLow, const double &pbcHigh) {
    // find the periodic image of pos with shortest distance to trg

    double pbc = pbcHigh - pbcLow;

    // find k s.t. min (pos-trg)^2 for pos = pos + k * pbc
    int k = int(round((trg - pos) / pbc));
    pos = pos + k * pbc;
}

void makeOutputDirectory(char *dir_name) {
    struct stat st;
    if (stat(dir_name, &st) != 0) {
        PS::S32 ret_loc = 0;
        PS::S32 ret = 0;
        if (PS::Comm::getRank() == 0)
            ret_loc = mkdir(dir_name, 0777);
        PS::Comm::broadcast(&ret_loc, ret);
        if (ret == 0) {
            if (PS::Comm::getRank() == 0)
                fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
        } else {
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
            PS::Abort();
        }
    }
}

void setBoundary(PS::DomainInfo &dinfo, const Config &runConfig) {
    if (runConfig.shell == true) {
        // use open boundary only
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    } else {
        const int pbcX = (runConfig.xPeriodic > 0 ? 1 : 0);
        const int pbcY = (runConfig.yPeriodic > 0 ? 1 : 0);
        const int pbcZ = (runConfig.zPeriodic > 0 ? 1 : 0);
        const int pbcFlag = 100 * pbcX + 10 * pbcY + pbcZ;

        switch (pbcFlag) {
            case 0:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
                break;
            case 1:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
                break;
            case 10:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
                break;
            case 100:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
                break;
            case 11:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
                break;
            case 101:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
                break;
            case 110:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
                break;
            case 111:
                dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
                break;
        }
    }
    dinfo.setPosRootDomain(PS::F64vec3(runConfig.simBoxLow[0], runConfig.simBoxLow[1], runConfig.simBoxLow[2]),
            PS::F64vec3(runConfig.simBoxHigh[0], runConfig.simBoxHigh[1], runConfig.simBoxHigh[2]));
    // rootdomain must be specified after PBC
}

// for the given change rate in unit time, calculate whether the change event happen in the current time step
inline bool calcChange(const double changeRateUnitTime, const double dt, const double rngU01) {
    // assume Poisson process.
    double pChangeDt = changeRateUnitTime * dt;
    return rngU01 < pChangeDt ? true : false;
}

bool readPosInitial(std::string posFile, std::vector<SphereFP> &sphereRead, int &nSphere) {
    std::ifstream myfile(posFile);
    if (!myfile.is_open()) {
        printf("Error: cannot open initial position file. use random initial configuration\n");
        sphereRead.clear();
        return false;
    }
    std::string line;
    getline(myfile, line); //
    getline(myfile, line); // read two header lines

    nSphere = 0;

    while (std::getline(myfile, line)) {
        std::istringstream liness(line);
        SphereFP newBody;
        nSphere++;
        int gid;
        char type;
        double px, py, pz;
        double radius;
        liness >> type >> gid >> radius >> px >> py >> pz;

        newBody.pos[0] = px;
        newBody.pos[1] = py;
        newBody.pos[2] = pz;
        newBody.gid = gid;
        newBody.radius = radius;
        std::cout << gid << " " << radius << " " << px << " " << py << " " << " " << pz << std::endl;
        sphereRead.push_back(newBody);
    }
    myfile.close();

    std::cout << "read initial finished" << std::endl;
    return true;
}

void SphereSystem::setInit(std::string posFile) {

    // RNG pool is not initialized yet
    // this function executes only on process 0
    // use c++ std rng instead of rng pool

    double minBoxEdge = std::min(runConfig.simBoxHigh[0] - runConfig.simBoxLow[0],
            runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]);
    minBoxEdge = std::min(runConfig.simBoxHigh[0] - runConfig.simBoxLow[0], minBoxEdge);

    // rng generator for polydispersed head and tail
    std::unique_ptr<TRngPool> swimmerHeadRngPoolPtr(new TRngPool(runConfig.rngSeed + 1));
    std::unique_ptr<TRngPool> swimmerTailRngPoolPtr(new TRngPool(runConfig.rngSeed + 2));

    if (runConfig.sphereRadiusSigmaA > 0) {
        swimmerHeadRngPoolPtr->setLogNormalParameters(runConfig.sphereRadiusA, runConfig.sphereRadiusSigmaA);
    }
    if (runConfig.sphereRadiusSigmaB > 0) {
        swimmerTailRngPoolPtr->setLogNormalParameters(runConfig.sphereRadiusB, runConfig.sphereRadiusSigmaB);
    }

    std::vector<SphereFP> sphereRead;
    sphereRead.reserve(runConfig.sphereNumber);
    int nSphereFile = 0;
    bool readFromFile = readPosInitial(posFile, sphereRead, nSphereFile);

    int nSphereGlobal = (nSphereFile == 0 ? runConfig.sphereNumber : nSphereFile);
    if (nSphereGlobal != runConfig.sphereNumber) {
        std::cout << "number of tubules found in initial position file" << std::endl;
        std::cout << "does NOT match the objects number in runconfig file" << std::endl;
        std::cout << "using the actual number of tubules found in position file" << std::endl;
    }

    if (isOdd(nSphereGlobal)) {
        std::cout << "Global sphere number must be an even number" << std::endl;
        exit(1);
    }

    systemSP.setNumberOfParticleLocal(nSphereGlobal);

    int repeatIndex = 0;

    PS::F64vec3 boxHigh(runConfig.simBoxHigh[0], runConfig.simBoxHigh[1], runConfig.simBoxHigh[2]);
    PS::F64vec3 boxLow(runConfig.simBoxLow[0], runConfig.simBoxLow[1], runConfig.simBoxLow[2]);
    const auto center = (boxHigh + boxLow) * 0.5;
    const double radius = std::min(std::min(boxHigh.x - center.x, boxHigh.y - center.y), boxHigh.z - center.z);

    // step 1 initialize spheres
    double radiusMax = 0;
    for (int i = 0; i < nSphereGlobal; i++) {
        bool swimmerHead = !isOdd(i); // swimmer head id = 0,2,4,6,8,...

        if (!readFromFile) {

            // not from file
            systemSP[i].gid = i;

            // generate within box and spherical shell
            bool accept = false;
            while (accept == false) {
                if (!swimmerHead) {
                    PS::F64vec3 direction(0, 0, 0);
                    sphereRngUVU01(myRngPoolPtr->getU01(), myRngPoolPtr->getU01(), direction);
                    systemSP[i].pos = systemSP[i - 1].pos + runConfig.springLength * direction;
                } else {
                    systemSP[i].pos = PS::F64vec3(
                            (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0]) * myRngPoolPtr->getU01()
                                    + runConfig.simBoxLow[0],
                            (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) * myRngPoolPtr->getU01()
                                    + runConfig.simBoxLow[1],
                            (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]) * myRngPoolPtr->getU01()
                                    + runConfig.simBoxLow[2]);
                }
                if (runConfig.shell == false
                        || (runConfig.shell == true && normVec3(systemSP[i].pos - center) < radius)) {
                    accept = true;
                }
            }

            if (swimmerHead) {
                if (runConfig.sphereRadiusSigmaA > 0) {
                    systemSP[i].radius = swimmerHeadRngPoolPtr->getLN();
                } else {
                    systemSP[i].radius = runConfig.sphereRadiusA;
                }
            } else {
                if (runConfig.sphereRadiusSigmaB > 0) {
                    systemSP[i].radius = swimmerTailRngPoolPtr->getLN();
                } else {
                    systemSP[i].radius = runConfig.sphereRadiusB;
                }
            }

        } else {
            // from file
            // system behavior is undefined if the data in the file doesn't make sense for pairs of spheres
            assert(i < nSphereFile);

            systemSP[i].gid = sphereRead[i].gid;
            systemSP[i].pos = sphereRead[i].pos;
            systemSP[i].radius = sphereRead[i].radius;
        }
        radiusMax = std::max(radiusMax, systemSP[i].radius);

        systemSP[i].posT0 = PS::F64vec3(0, 0, 0);
        zeroVec3(systemSP[i].forceExt); // external...
        zeroVec3(systemSP[i].torqueExt);
        zeroVec3(systemSP[i].forceNear);
        zeroVec3(systemSP[i].torqueNear);
        zeroVec3(systemSP[i].forceCollision);
        zeroVec3(systemSP[i].torqueCollision);
        zeroVec3(systemSP[i].vel);
        zeroVec3(systemSP[i].omega);

        if (runConfig.monolayer == true) {
            const double monoZ = (runConfig.simBoxHigh[2] + runConfig.simBoxLow[2]) * 0.5;
            systemSP[i].pos[2] = monoZ;
        }
    }

    // set search radius
    std::cout << "max Radius: " << radiusMax << std::endl;
    for (int i = 0; i < systemSP.getNumberOfParticleLocal(); i++) {
        // for tubule
        systemSP[i].RSearch = radiusMax * (2 + 0.5);
    }

    // output for debug
    for (int i = 0; i < nSphereGlobal / 2; i++) {
        auto &head = systemSP[2 * i];
        auto &tail = systemSP[2 * i + 1];

        auto direction = head.pos - tail.pos;
        auto length = normVec3(direction);
        direction *= (1 / length);
        printf("swimmer -------------------------------------------\n");
        printf("id: %d, head radius: %.6lf, tail radius: %.6lf\n", i, head.radius, tail.radius);
        printf("direction: %lf,%lf,%lf\n", direction.x, direction.y, direction.z);
        printf("length: %lf, neighbor search radius: %lf\n", length, head.RSearch);
        printf("---------------------------------------------------\n");
    }

    // check volume fraction
    double particleVolume = 0;
    for (int i = 0; i < systemSP.getNumberOfParticleLocal(); i++) {
        const auto &sphere = systemSP[i];
        particleVolume += 3.1416 * (4 / 3.0) * pow(sphere.radius, 3);
    }
    double boxVolume = (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0])
            * (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) * (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]);

    std::cout << "initial volume fraction: " << particleVolume / boxVolume << std::endl;
}

SphereSystem::SphereSystem(std::string configFile, std::string posFile, int argc, char **argv) :
        runConfig(configFile), calcNearForceFtr() {
    // TRNG pool must be initialized after mpi is initialized

    const int pbcX = (runConfig.xPeriodic > 0 ? 1 : 0);
    const int pbcY = (runConfig.yPeriodic > 0 ? 1 : 0);
    const int pbcZ = (runConfig.zPeriodic > 0 ? 1 : 0);
    const int pbcFlag = 100 * pbcX + 10 * pbcY + pbcZ;

    // initialize system
    MPI_Barrier(MPI_COMM_WORLD);

    this->myRank = PS::Comm::getRank();
    this->nProcs = PS::Comm::getNumberOfProc();
    this->sysTime = 0;
    this->snapTime = 0;
    this->id_snap = 0;

    this->systemSP.initialize();

    // initialize rng pool after mpi is initialized
    // no make_unique in c++11 (introduced in c++14)
    // myRngPoolPtr = std::make_unique<TRngPool>(runConfig.rngSeed);
    myRngPoolPtr = std::move(std::unique_ptr<TRngPool>(new TRngPool(runConfig.rngSeed)));

    if (PS::Comm::getRank() == 0) {
        // prepare the output directory
        sprintf(this->dir_name, "./result");
        makeOutputDirectory(dir_name);
        // initialize on head rank
        setInit(posFile);
    } else {
        systemSP.setNumberOfParticleLocal(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // set domain info
    const PS::F64 coef_ema = 0.5;
    this->dinfo.initialize(coef_ema);
    setBoundary(this->dinfo, this->runConfig);
    this->systemSP.adjustPositionIntoRootDomain(dinfo);
    this->dinfo.decomposeDomainAll(systemSP);

    // initial exchange of particles
    MPI_Barrier(MPI_COMM_WORLD);
    systemSP.exchangeParticle(dinfo);

    printf("n body local: %d on process %d\n", systemSP.getNumberOfParticleLocal(), PS::Comm::getRank());

    printf("initial configuration initialized, process %d\n", PS::Comm::getRank());

    // initialize force calculate tree
    const PS::U64 ntot = systemSP.getNumberOfParticleGlobal();
    std::cout << "ntot: " << ntot << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // be careful, strange bugs may appear if tree is not initialized with large enough n_glob
    this->treeNear.initialize(ntot * 4);

    PS::F64vec3 boxHigh(runConfig.simBoxHigh[0], runConfig.simBoxHigh[1], runConfig.simBoxHigh[2]);
    PS::F64vec3 boxLow(runConfig.simBoxLow[0], runConfig.simBoxLow[1], runConfig.simBoxLow[2]);

    if (runConfig.shell == true) {
        auto center = (boxHigh + boxLow) * 0.5;
        const double radius = std::min(std::min(boxHigh.x - center.x, boxHigh.y - center.y), boxHigh.z - center.z);
        shell.initialize(radius, center);
    }

    // initialize hydrosystem

    if (runConfig.hydro == true) {
        const int mp_order = runConfig.pFMM;
        const double mu = runConfig.viscosity;
        std::array<double, 3> boxHigh = { runConfig.simBoxHigh[0], runConfig.simBoxHigh[1], runConfig.simBoxHigh[2] };
        std::array<double, 3> boxLow = { runConfig.simBoxLow[0], runConfig.simBoxLow[1], runConfig.simBoxLow[2] };
        FMM_Wrapper::PAXIS pbc;
        switch (pbcFlag) {
            case 0:
                std::cout << "using free boundary non-periodic box" << std::endl;
                pbc = FMM_Wrapper::PAXIS::NONE;
                break;
            case 1:
                std::cout << "using 1D periodic box PZ" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PZ;
                break;
            case 10:
                std::cout << "using 1D periodic box PY" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PY;
                break;
            case 100:
                std::cout << "using 1D periodic box PX" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PX;
                break;
            case 11:
                std::cout << "using 2D periodic box PYZ" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PYZ;
                break;
            case 101:
                std::cout << "using 2D periodic box PXZ" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PXZ;
                break;
            case 110:
                std::cout << "using 2D periodic box PXY" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PXY;
                break;
            case 111:
                std::cout << "using 3D periodic box" << std::endl;
                pbc = FMM_Wrapper::PAXIS::PXYZ;
                break;
        }
        pointSphereSystemPtr.reset(new HydroSphereSystem(mu, pbc, mp_order, boxLow, boxHigh));
    }

    // initial resolve and reset move to zero
    // 5 pass to resolve initial overlaps
    this->withHydro = false;
    double bufferGap = systemSP[0].RSearch * 0.2;
    for (int i = 0; i < 5; i++) {
        prepareT0();
        treeNear.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);
        calcSpringForce();
        moveAll(bufferGap);
        bufferGap *= 0.5;
    }
    output();

    this->withHydro = runConfig.hydro;

    printf("SphereSystem initialized on process: %d \n", PS::Comm::getRank());
}

// first order Euler
void SphereSystem::moveAll(double bufferGap) {
    const double dt = runConfig.dt; // half step for mid-point
    const double Pi = 3.14159265358979323846;
    const double mu = runConfig.viscosity;
    const double scaleBrown = runConfig.scaleBrown;
    const double kBT = runConfig.kBT;

    const int nLocal = systemSP.getNumberOfParticleLocal();
    const auto omega = PS::F64vec3(0, 1, 0);
    PS::F64vec3 boxHigh(runConfig.simBoxHigh[0], runConfig.simBoxHigh[1], runConfig.simBoxHigh[2]);
    PS::F64vec3 boxLow(runConfig.simBoxLow[0], runConfig.simBoxLow[1], runConfig.simBoxLow[2]);
    const auto center = (boxHigh + boxLow) * 0.5;
    const auto dist = boxHigh - center;
    const double radius = std::min(std::min(boxHigh.x - center.x, boxHigh.y - center.y), boxHigh.z - center.z);

// step 1 : calc velocity from external, Brownian and pair force
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        SphereFP &sphere = systemSP[i];
        const double radius = sphere.radius;
        const double drag = 6 * Pi * radius * mu;
        const double dragRot = 8 * Pi * mu * radius * radius * radius;

        const auto FBtrans = scaleBrown * sphere.randTrans * sqrt(2 * kBT * drag / (dt));
        const auto TBrot = scaleBrown * sphere.randRot * sqrt(2 * kBT * dragRot / (dt));

        sphere.randTrans = FBtrans;
        sphere.randRot = TBrot;

        const auto forceTotal = sphere.forceExt + FBtrans + sphere.forceNear;
        const auto torqueTotal = sphere.torqueExt + TBrot + sphere.torqueNear;

        sphere.vel = forceTotal / drag;
        sphere.omega = torqueTotal / dragRot;

#ifdef MYDEBUGINFO
        std::cout << "force total" << forceTotal << std::endl;
        std::cout << "torque total" << torqueTotal << std::endl;
        std::cout << "vel: " << sphere.vel << std::endl;
        std::cout << "omega: " << sphere.omega << std::endl;
#endif
    }

    // step 1.2: calc velocity from hydro.
    // Hydro is calculated with Fnear. Fcol is implicit
    if (this->withHydro == true && runConfig.hydro == true) {
        pointSphereSystemPtr->pointSphere.resize(nLocal);

#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            auto &pointSphere = pointSphereSystemPtr->pointSphere[i];
            const auto &sphere = systemSP[i];
            for (int j = 0; j < 3; j++) {
                pointSphere.pos[j] = sphere.pos[j];
                pointSphere.radius = sphere.radius;
                pointSphere.force[j] = sphere.forceNear[j];
                // key assumption of model:
                // forceExt = fswim does not go into FMM
                //                std::cout << "point force" << pointSphere.force[j] << std::endl;
            }
        }

        pointSphereSystemPtr->updateSystem();
        pointSphereSystemPtr->farfieldMobilityApply();

#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            const auto &pointSphere = pointSphereSystemPtr->pointSphere[i];
            auto &sphere = systemSP[i];
            for (int j = 0; j < 3; j++) {
                sphere.vel[j] += pointSphere.vel[j];
            }
        }
    }

    // step 2 resolve collision force with LCP
    // for implicit stepping, vel due to collision is calculated in this step
    // overlapResolve(bufferGap);

// step 3 apply move;
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        SphereFP &sphere = systemSP[i];

        sphere.pos += sphere.vel * dt;

#ifdef MYDEBUGINFO
        std::cout << "vel: " << sphere.vel << std::endl;
        std::cout << "omega: " << sphere.omega << std::endl;
#endif
    }

    // save force to hydrosystem to generate flow
    if (this->withHydro == true && runConfig.hydro == true) {
        pointSphereSystemPtr->pointSphere.resize(nLocal);

#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            auto &pointSphere = pointSphereSystemPtr->pointSphere[i];
            const auto &sphere = systemSP[i];
            for (int j = 0; j < 3; j++) {
                pointSphere.pos[j] = sphere.pos[j];
                pointSphere.radius = sphere.radius;
                pointSphere.force[j] = sphere.forceNear[j] + sphere.forceCollision[j];
                // key assumption of model:
                // forceExt = fswim does not go into FMM
                //                std::cout << "point force" << pointSphere.force[j] << std::endl;
            }
        }
    }
    if (runConfig.monolayer == true) {
        const double monoZ = 0.5 * (runConfig.simBoxHigh[2] + runConfig.simBoxLow[2]);
#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            systemSP[i].pos.z = monoZ;
        }
    }
}

void SphereSystem::stepForward() {
    // determine output
    prepareT0();

    if ((sysTime >= snapTime) || ((sysTime + runConfig.dt) - snapTime) > (snapTime - sysTime)) {
        dumpxyz = true;
        dumpflow = true;
        id_snap++;
        snapTime += runConfig.timeSnap;
    }

    output(); // guarantee the structure and flow field are generated with the same configuration.

// step 1: generate RNGs for each object
// for tubule, generate 7N01 for rotation and translation
// for protein, generate 2 U01 for bind, 2 U01 for orientation, 3 N01 for translation
#ifdef MYDEBUGINFO
    printf("T0 prepared: %d \n", PS::Comm::getRank());
    std::ofstream myfile;
    std::stringstream rankNumber;
    rankNumber << "Treedump_" << PS::Comm::getRank();
    myfile.open(rankNumber.str().c_str());
    this->treeNear.dump(myfile);
    myfile.close();
#endif

    // step 2: Euler move
    treeNear.calcForceAllAndWriteBack(calcNearForceFtr, systemSP, dinfo);
    calcSpringForce();
    moveAll(0);

#ifdef MYDEBUGINFO
    printf("full step moved: %d \n", PS::Comm::getRank());
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    statistics();

    systemSP.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(systemSP);

    MPI_Barrier(MPI_COMM_WORLD);
    systemSP.exchangeParticle(dinfo);

    printf("%lf\n", sysTime);
    sysTime += runConfig.dt;

#ifdef MYDEBUGINFO
    printf("particle exchanged: %d \n", PS::Comm::getRank());
#endif
}

void SphereSystem::output() {
    //	if ((sysTime >= snapTime) || ((sysTime + runConfig.dt) - snapTime) > (snapTime - sysTime)) {
    if (dumpxyz) {
        dumpxyz = false;
        char filename[256];
        sprintf(filename, "%s/%05d.xyz", dir_name, id_snap);
        FileHeader header;
        header.time = sysTime;
        header.nSphere = systemSP.getNumberOfParticleGlobal();
        MPI_Barrier(MPI_COMM_WORLD);
        systemSP.writeParticleAscii(filename, header);

        if (pointSphereSystemPtr && runConfig.hydro && runConfig.dumpflow) {
            char flowfile[256];
            sprintf(flowfile, "%s/FLOW%05d", dir_name, id_snap);
            std::string flowFileName(flowfile);
            pointSphereSystemPtr->calcFlowOnGrid(0.1, flowFileName, 0, 0, 0);
        }
    }
}

void SphereSystem::prepareT0() {
    const double Pi = 3.14159265358979323846;
    const double mu = runConfig.viscosity;
    const double dt = runConfig.dt;

    // clear nearforcecalculator
    for (auto &colBlockQue : *((this->calcNearForceFtr).colBlockThread_ptr)) {
        colBlockQue.clear();
    }

    // create a map, get global indices
    const int nSphereLocal = systemSP.getNumberOfParticleLocal();
    const int nSphereGlobal = systemSP.getNumberOfParticleGlobal();
    Teuchos::RCP<const TCOMM> commRcp = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    this->myRank = commRcp->getRank();
    this->nProcs = commRcp->getSize();
    // set up contiguous map with global and local number
    Teuchos::RCP<const TMAP> spMapRcp = Teuchos::rcp(new TMAP(nSphereGlobal, nSphereLocal, 0, commRcp));
    const int minGlobalIndex = spMapRcp->getMinGlobalIndex();
    // TRNGpool is thread safe
    // supply known threadId for better speed
    // generate RNG and save postition T0

    const int nTotal = this->systemSP.getNumberOfParticleLocal();
    assert(nSphereLocal == nTotal);
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (size_t i = 0; i < nTotal; i++) {
            auto &sphere = systemSP[i];
            sphere.overlap = 0;
            sphere.localIndex = i;
            sphere.globalIndex = minGlobalIndex + i;
            sphere.randTrans[0] = myRngPoolPtr->getN01(threadId); // translational Brownian force
            sphere.randTrans[1] = myRngPoolPtr->getN01(threadId);
            sphere.randTrans[2] = myRngPoolPtr->getN01(threadId);
            sphere.randRot[0] = myRngPoolPtr->getN01(threadId);
            sphere.randRot[1] = myRngPoolPtr->getN01(threadId);
            sphere.randRot[2] = myRngPoolPtr->getN01(threadId); // this component of RNG is not used
            zeroVec3(sphere.forceExt); // set zero initial
            zeroVec3(sphere.torqueExt); // set zero initial
            zeroVec3(sphere.forceNear); // set zero initial
            zeroVec3(sphere.torqueNear); // set zero initial
            zeroVec3(sphere.vel); // set zero initial
            zeroVec3(sphere.omega); // set zero initial
            zeroVec3(sphere.forceCollision);
            zeroVec3(sphere.torqueCollision);
            // save T0
            sphere.posT0 = sphere.pos;
            // forceCollision should be preserved from last timestep to use explicit timestepping
        }
    }
}

/*
 Resolve overlap with LCP solver
 Precondition: 1. collision pairs colBlockThread_ptr,
 2. non-collision velocities in sphere.vel
 Postcondition: 1. collision force magnitude(lambda) in colBlockThread_ptr for stress calculation
 2. collision force in sphere.forceCollision
 Note: for spheres, omegas/torques do not appear in collision resolving
 */
void SphereSystem::overlapResolve(double bufferGap) {
    const double tol = 1e-5;

    if (runConfig.shell == true) {
        shell.clear();
        const int ntotal = systemSP.getNumberOfParticleLocal();
#pragma omp parallel
        {
            const int threadId = omp_get_thread_num();
#pragma omp for
            for (int is = 0; is < ntotal; is++) {
                auto &sphere = systemSP[is];
                shell.processParticle(sphere, threadId);
            }
        }
    } else {
        // clear boundary colblocks
        shell.clear();
    }

#ifdef DEBUGLCPCOL
    // dump collision data
    for (auto &colBlockQue : *((this->calcNearForceFtr).colBlockThread_ptr)) {
        std::cout << colBlockQue.size() << "particle collisions" << std::endl;
        for (auto &colBlock : colBlockQue) {
            // std::cout << colBlock.gidI << " " << colBlock.gidJ << "  sep:" << colBlock.phi0 << std::endl;
            std::cout << colBlock.globalIndexI << " " << colBlock.globalIndexJ << "  sep:" << colBlock.phi0
            << std::endl;
        }
    }

    // dump collision data
    for (auto &colBlockQue : *(shell.colBlockThread_ptr)) {
        std::cout << colBlockQue.size() << "shell collisions" << std::endl;
        for (auto &colBlock : colBlockQue) {
            // std::cout << colBlock.gidI << " " << colBlock.gidJ << "  sep:" << colBlock.phi0 << std::endl;
            std::cout << colBlock.globalIndexI << " " << colBlock.globalIndexJ << "  sep:" << colBlock.phi0
            << std::endl;
        }
    }

#endif

    Teuchos::RCP<const TCOMM> commRcp = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    // step 1 create maps
    const int nSphereLocal = systemSP.getNumberOfParticleLocal();
    const int nSphereGlobal = systemSP.getNumberOfParticleGlobal();
    // sphere map contiguous map with global and local number
    Teuchos::RCP<const TMAP> sphereMapRcp = Teuchos::rcp(new TMAP(nSphereGlobal, nSphereLocal, 0, commRcp));
    Teuchos::RCP<const TMAP> mobMatMapRcp = Teuchos::rcp(new TMAP(3 * nSphereGlobal, 3 * nSphereLocal, 0, commRcp));

    // collision pair map
    // particle-particle collisions
    const auto &colBlockThread = *(this->calcNearForceFtr.colBlockThread_ptr);
    const int nThreads = colBlockThread.size();
    std::vector<int> colBlockThreadIndex(nThreads + 1);
    colBlockThreadIndex[0] = 0;
    for (int i = 0; i < nThreads; i++) {
        colBlockThreadIndex[i + 1] = colBlockThreadIndex[i] + colBlockThread[i].size();
    }
    // particle-shell collisions
    const auto &shellBlockThread = *(shell.colBlockThread_ptr);
    std::vector<int> shellBlockThreadIndex(nThreads + 1);
    shellBlockThreadIndex[0] = 0;
    for (int i = 0; i < nThreads; i++) {
        shellBlockThreadIndex[i + 1] = shellBlockThreadIndex[i] + shellBlockThread[i].size();
    }
    // global reduce
    const int localColBlock = colBlockThreadIndex.back() + shellBlockThreadIndex.back();
    int globalColBlock = localColBlock;
    MPI_Allreduce(MPI_IN_PLACE, &globalColBlock, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::cout << " globalColBlock " << globalColBlock << std::endl;
    std::cout << " localColBlock " << localColBlock << std::endl;
    Teuchos::RCP<const TMAP> colBlockMapRcp = Teuchos::rcp(new TMAP(globalColBlock, localColBlock, 0, commRcp));

    // phi0 and Vn vectors
    Teuchos::RCP<TV> phi0VecRcp = Teuchos::rcp(new TV(colBlockMapRcp.getConst(), false));
    Teuchos::RCP<TV> vnVecRcp = Teuchos::rcp(new TV(mobMatMapRcp.getConst(), false));
    // raw ptrs to local values
    auto phi0Vec_2d = phi0VecRcp->getLocalView<Kokkos::HostSpace>();    // LeftLayout
    phi0VecRcp->modify<Kokkos::HostSpace>();
    const double dt = runConfig.dt;

    const double maxVel = 10.0; // a 'stabilizer', limit the maximum moving apart speed to be maxVel um/s, avoid unbounded collision force
    // ref. Page 24 of Tasora, A. Time integration in ChronoEngine. (2016).
    // not used here.

#pragma omp parallel for
    for (int threadId = 0; threadId < nThreads; threadId++) {
        // particle-particle
        const auto &colBlockMyThread = colBlockThread[threadId];
        const int colBlockNum = colBlockMyThread.size();
        const int colBlockIndexBase = colBlockThreadIndex[threadId];
        for (int j = 0; j < colBlockNum; j++) {
            const double phi0buf = colBlockMyThread[j].phi0 / dt - bufferGap;
            phi0Vec_2d(colBlockIndexBase + j, 0) = phi0buf;    // std::max(phi0buf, -maxVel);
        }
        // particle-boundary
        const auto &shellBlockMyThread = shellBlockThread[threadId];
        const int shellBlockNum = shellBlockMyThread.size();
        const int shellBlockIndexBase = shellBlockThreadIndex[threadId] + colBlockThreadIndex.back();
        for (int j = 0; j < shellBlockNum; j++) {
            const double phi0buf = shellBlockMyThread[j].phi0 / dt - bufferGap;
            phi0Vec_2d(shellBlockIndexBase + j, 0) = phi0buf;    // std::max(phi0buf, -maxVel);
        }
    }

    std::cout << " phi0vec filled " << std::endl;

    // initialize the Fc^T collision matrix
    // rowMap = colBlockMap, colMap=mobMatMap
    // six entries for each row (each colBlock), gI[x,y,z], gJ[x,y,z]
    Kokkos::View<size_t *> rowPointers("rowPointers", localColBlock + 1);
    rowPointers[0] = 0;
    // particle-particle blocks
    for (int i = 1; i <= colBlockThreadIndex.back(); i++) {
        rowPointers[i] = rowPointers[i - 1] + 6;
    }
    // particle-shell blocks
    for (int i = colBlockThreadIndex.back() + 1; i <= shellBlockThreadIndex.back() + colBlockThreadIndex.back(); i++) {
        rowPointers[i] = rowPointers[i - 1] + 3;
    }
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localColBlock]);
    Kokkos::View<double *> values("values", rowPointers[localColBlock]);
#pragma omp parallel for num_threads(nThreads)
    // get the column indices
    for (int threadId = 0; threadId < nThreads; threadId++) {
        // particle-particle collision
        const auto &colBlockMyThread = colBlockThread[threadId];
        const int colBlockNum = colBlockMyThread.size();
        const int colBlockIndexBase = colBlockThreadIndex[threadId];
        for (int j = 0; j < colBlockNum; j++) {
            const int kk = 6 * (colBlockIndexBase + j);
            columnIndices[kk] = 3 * colBlockMyThread[j].globalIndexI;
            columnIndices[kk + 1] = 3 * colBlockMyThread[j].globalIndexI + 1;
            columnIndices[kk + 2] = 3 * colBlockMyThread[j].globalIndexI + 2;
            columnIndices[kk + 3] = 3 * colBlockMyThread[j].globalIndexJ;
            columnIndices[kk + 4] = 3 * colBlockMyThread[j].globalIndexJ + 1;
            columnIndices[kk + 5] = 3 * colBlockMyThread[j].globalIndexJ + 2;
            values[kk] = colBlockMyThread[j].gvecI[0];
            values[kk + 1] = colBlockMyThread[j].gvecI[1];
            values[kk + 2] = colBlockMyThread[j].gvecI[2];
            values[kk + 3] = colBlockMyThread[j].gvecJ[0];
            values[kk + 4] = colBlockMyThread[j].gvecJ[1];
            values[kk + 5] = colBlockMyThread[j].gvecJ[2];
        }
        // particle-shell collision
        const auto &shellBlockMyThread = shellBlockThread[threadId];
        const int shellBlockNum = shellBlockMyThread.size();
        const int shellBlockIndexBase = shellBlockThreadIndex[threadId];
        for (int j = 0; j < shellBlockNum; j++) {
            const int kk = 3 * (shellBlockIndexBase + j) + 6 * colBlockThreadIndex.back();
            // std::cout << kk << std::endl;
            // std::cout << j << std::endl;
            // std::cout << shellBlockIndexBase << std::endl;
            // std::cout << shellBlockNum << std::endl;
            columnIndices[kk] = 3 * shellBlockMyThread[j].globalIndexI;
            columnIndices[kk + 1] = 3 * shellBlockMyThread[j].globalIndexI + 1;
            columnIndices[kk + 2] = 3 * shellBlockMyThread[j].globalIndexI + 2;
            values[kk] = shellBlockMyThread[j].gvecI[0];
            values[kk + 1] = shellBlockMyThread[j].gvecI[1];
            values[kk + 2] = shellBlockMyThread[j].gvecI[2];
        }
    }

    std::vector<int> globalIndexOnLocal(3 * nSphereGlobal);
#pragma omp parallel for
    // define the matrix with col map = full cols on every node
    for (int kk = 0; kk < 3 * nSphereGlobal; kk++) {
        globalIndexOnLocal[kk] = kk;
    }
    Teuchos::RCP<TMAP> colMapRcp = Teuchos::rcp(
            new TMAP(3 * nSphereGlobal, globalIndexOnLocal.data(), 3 * nSphereGlobal, 0, commRcp));

    Teuchos::RCP<TCMAT> fcTransMatRcp = Teuchos::rcp(
            new TCMAT(colBlockMapRcp, colMapRcp, rowPointers, columnIndices, values));
    fcTransMatRcp->fillComplete(mobMatMapRcp, colBlockMapRcp); // domainMap, rangeMap

#ifdef DEBUGLCPCOL
    std::cout << "FcTransConstructed: " << fcTransMatRcp->description() << std::endl;
    dumpTCMAT(fcTransMatRcp, "fcTransMat.mtx");
#endif

    auto vn_2d = vnVecRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    vnVecRcp->modify<Kokkos::HostSpace>();
    for (int i = 0; i < nSphereLocal; i++) {
        vn_2d(3 * i, 0) = systemSP[i].vel[0];
        vn_2d(3 * i + 1, 0) = systemSP[i].vel[1];
        vn_2d(3 * i + 2, 0) = systemSP[i].vel[2];
    }

    // the vec b in LCP stored in phi0
    // phi0 = phi0 + Fc^T * Vn
    MPI_Barrier(MPI_COMM_WORLD);
    fcTransMatRcp->apply(*vnVecRcp, *phi0VecRcp, Teuchos::NO_TRANS, 1.0, 1.0);
#ifdef DEBUGLCPCOL
    dumpTMV(phi0VecRcp, "phi0VecRcp.mtx");
#endif
    int phi0NegNumber = 0;
#pragma omp parallel for reduction(+ : phi0NegNumber)
    for (int j = 0; j < phi0Vec_2d.dimension_0(); j++) {
        if (phi0Vec_2d(j, 0) < -tol) {
            phi0NegNumber = phi0NegNumber + 1;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &phi0NegNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (myRank == 0) {
        std::cout << "neg phi0: " << phi0NegNumber << std::endl;
    }
    if (phi0NegNumber <= 0) {
        return;
    }

    Teuchos::RCP<TOP> mobOpRcp;

    if (runConfig.hydro == false || withHydro == false) {
        // mobility matrix operator
        const double Pi = 3.14159265358979323846;
        const double mu = runConfig.viscosity;
        // diagonal hydro mobility operator
        const int localSize = nSphereLocal * 3;
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
        for (int i = 0; i < nSphereLocal; i++) {
            const auto &sphere = systemSP[i];
            const double radius = sphere.radius;
            const double dragInv = 1.0 / (6 * Pi * radius * mu);
            values[3 * i] = dragInv;
            values[3 * i + 1] = dragInv;
            values[3 * i + 2] = dragInv;
        }
        Teuchos::RCP<TCMAT> mobMatRcp = Teuchos::rcp(
                new TCMAT(mobMatMapRcp, mobMatMapRcp, rowPointers, columnIndices, values));
        mobMatRcp->fillComplete(mobMatMapRcp, mobMatMapRcp); // domainMap, rangeMap
#ifdef DEBUGLCPCOL
                std::cout << "MobMatConstructed: " << mobMatRcp->description() << std::endl;
                dumpTCMAT(mobMatRcp, "MobMat.mtx");
#endif
        mobOpRcp = mobMatRcp;
    } else {
        // mobility matrix operator with full hydro
        // the FMM tree should have been udpated in moveAll()
        this->pointSphereSystemPtr->buildMobilityOperator();
        mobOpRcp = this->pointSphereSystemPtr->mobilityOperator;

    }
    // the operator A in LCP
    Teuchos::RCP<CPMatOp> AmatRcp = Teuchos::rcp(new CPMatOp(mobOpRcp, fcTransMatRcp));

    // create the solver
    IteHistory history;
    Teuchos::RCP<TV> xsolRcp = Teuchos::rcp(new TV(colBlockMapRcp.getConst(), true));        // zero initial guess
    const int maxIte = 2000;

    CPSolver myLCPSolver(AmatRcp, phi0VecRcp, fcTransMatRcp->getRangeMap(), commRcp);

    myLCPSolver.LCP_BBPGD(xsolRcp, tol, maxIte, history);

    if (myRank == 0) {
        auto &p = history.back();
        std::cout << "LCP residue: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << std::endl;
    }
    Teuchos::RCP<TV> forceColRcp = AmatRcp->forceVecRcp;
#ifdef DEBUGLCPCOL
    if (myRank == 0) {
        for (auto &p : history) {
            std::cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << std::endl;
        }
    }
    dumpTMV(forceColRcp, "forceCol.mtx");
    dumpTMV(xsolRcp, "xsol.mtx");
#endif

    // get the result collision force and torque
    const auto forceCol_2d = forceColRcp->getLocalView<Kokkos::HostSpace>();        // LeftLayout
    const auto velCol_2d = AmatRcp->velVecRcp->getLocalView<Kokkos::HostSpace>();

    // put a force limiter
    const double forceLimiter = runConfig.swimForce * runConfig.forceMaxRatio;
// store results for force/displacement calculation
#pragma omp parallel for
    for (int i = 0; i < nSphereLocal; i++) {
        auto & forceCollision = systemSP[i].forceCollision;
        const double fx = forceCol_2d(3 * i, 0);
        const double fy = forceCol_2d(3 * i + 1, 0);
        const double fz = forceCol_2d(3 * i + 2, 0);
        forceCollision.x = std::copysign(std::min(std::abs(fx), forceLimiter), fx);
        forceCollision.y = std::copysign(std::min(std::abs(fy), forceLimiter), fy);
        forceCollision.z = std::copysign(std::min(std::abs(fz), forceLimiter), fz);
        auto & vel = systemSP[i].vel;
        vel.x += velCol_2d(3 * i, 0);
        vel.y += velCol_2d(3 * i + 1, 0);
        vel.z += velCol_2d(3 * i + 2, 0);
//        systemSP[i].forceCollision[0] = forceCol_2d(3 * i, 0);
//        systemSP[i].forceCollision[1] = forceCol_2d(3 * i + 1, 0);
//        systemSP[i].forceCollision[2] = forceCol_2d(3 * i + 2, 0);
    }

    // store results for stress calculation
    const auto xsol_2d = xsolRcp->getLocalView<Kokkos::HostSpace>();        // LeftLayout
#pragma omp parallel for num_threads(nThreads)
    for (int threadId = 0; threadId < nThreads; threadId++) {
        auto &colBlockMyThread = (*(this->calcNearForceFtr.colBlockThread_ptr))[threadId];
        const int colBlockNum = colBlockMyThread.size();
        const int colBlockIndexBase = colBlockThreadIndex[threadId];
        for (int j = 0; j < colBlockNum; j++) {
            colBlockMyThread[j].lambda = xsol_2d(colBlockIndexBase + j, 0);
        }
    }
    return;
}

void SphereSystem::statistics() {
    const double dt = runConfig.dt; // half step for mid-point
    const double Pi = 3.14159265358979323846;
    const double mu = runConfig.viscosity;
    const double scaleBrown = runConfig.scaleBrown;
    const double kBT = runConfig.kBT;

    const int nSphereLocal = systemSP.getNumberOfParticleLocal();
    const int nSphereGlobal = systemSP.getNumberOfParticleGlobal();
    const int nThreads = (this->calcNearForceFtr.colBlockThread_ptr)->size();

    double xx = 0, xy = 0, xz = 0, yx = 0, yy = 0, yz = 0, zx = 0, zy = 0, zz = 0;
// accumulate XF
#pragma omp parallel for num_threads(nThreads), reduction(+ : xx, xy, xz, yx, yy, yz, zx, zy, zz)
    for (int threadId = 0; threadId < nThreads; threadId++) {
        auto &colBlockMyThread = (*(this->calcNearForceFtr.colBlockThread_ptr))[threadId];
        const int colBlockNum = colBlockMyThread.size();
        for (int j = 0; j < colBlockNum; j++) {
            const PS::F64vec3 f = colBlockMyThread[j].lambda * colBlockMyThread[j].gvecI;
            const PS::F64vec3 x = colBlockMyThread[j].posI - colBlockMyThread[j].posJ;
            xx = xx + x[0] * f[0];
            xy = xy + x[0] * f[1];
            xz = xz + x[0] * f[2];
            yx = yx + x[1] * f[0];
            yy = yy + x[1] * f[1];
            yz = yz + x[1] * f[2];
            zx = zx + x[2] * f[0];
            zy = zy + x[2] * f[1];
            zz = zz + x[2] * f[2];
        }
    }

    // accumulate XF on all ranks
    double collisionXFLocal[9] = { xx, xy, xz, yx, yy, yz, zx, zy, zz };
    double collisionXF[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // MPI Reduce
    MPI_Allreduce(collisionXFLocal, collisionXF, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const double boxVolume = (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0])
            * (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) * (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]);
    const double numberDensity = systemSP.getNumberOfParticleGlobal() / boxVolume;
    const double nkBT = numberDensity * kBT;

    if (myRank == 0) {
        printf("collision Stress /(nKBT): %.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n",
                collisionXF[0] / (nSphereGlobal * kBT), collisionXF[1] / (nSphereGlobal * kBT),
                collisionXF[2] / (nSphereGlobal * kBT), collisionXF[3] / (nSphereGlobal * kBT),
                collisionXF[4] / (nSphereGlobal * kBT), collisionXF[5] / (nSphereGlobal * kBT),
                collisionXF[6] / (nSphereGlobal * kBT), collisionXF[7] / (nSphereGlobal * kBT),
                collisionXF[8] / (nSphereGlobal * kBT));
    }
}

void SphereSystem::calcSpringForce() {

    // calculate the spring force between head and tail.
    // for each head, find the tail position
    // for each tail, find the head position
    const int parNumber = systemSP.getNumberOfParticleLocal();
    ZDD<PS::F64vec3> dataFinder(parNumber);
    // allocate
    dataFinder.localID.resize(parNumber);
    dataFinder.localData.resize(parNumber);
// fill in data
#pragma omp parallel for schedule(dynamic, 512)
    for (int i = 0; i < parNumber; i++) {
        const auto &sphere = systemSP[i];
        dataFinder.localID[i] = sphere.gid;
        dataFinder.localData[i] = sphere.pos;
    }
    dataFinder.buildIndex();

    // specify search target
    // build id to find
    dataFinder.findID.resize(parNumber);
    dataFinder.findData.resize(parNumber);
#pragma omp parallel for schedule(dynamic, 512)
    for (int i = 0; i < parNumber; i++) {
        const auto &sphere = systemSP[i];
        bool swimmerHead = !isOdd(sphere.gid);
        dataFinder.findID[i] = swimmerHead ? sphere.gid + 1 : sphere.gid - 1; // head id 0 look for 1, tail id 1 look for 0, etc
    }

    MPI_Barrier(MPI_COMM_WORLD);
    dataFinder.find();

    if (runConfig.xPeriodic == 1) {
#pragma omp parallel for schedule(dynamic, 512)
        for (int i = 0; i < parNumber; i++) {
            periodicWrap(dataFinder.findData[i].x, systemSP[i].pos.x, runConfig.simBoxLow[0], runConfig.simBoxHigh[0]);
        }
    }

    if (runConfig.yPeriodic == 1) {
#pragma omp parallel for schedule(dynamic, 512)
        for (int i = 0; i < parNumber; i++) {
            periodicWrap(dataFinder.findData[i].y, systemSP[i].pos.y, runConfig.simBoxLow[1], runConfig.simBoxHigh[1]);
        }
    }

    if (runConfig.zPeriodic == 1) {
#pragma omp parallel for schedule(dynamic, 512)
        for (int i = 0; i < parNumber; i++) {
            periodicWrap(dataFinder.findData[i].z, systemSP[i].pos.z, runConfig.simBoxLow[2], runConfig.simBoxHigh[2]);
        }
    }

    // calc spring force
    // same notation as in the paper
    const auto &l = runConfig.springLength;
    const auto &lm = runConfig.springMaxLength;
    const double ratio = lm / l;
    const double Fmax = runConfig.swimForce * 20;

    int lmN = 0;
#pragma omp parallel for reduction(+:lmN)
    for (int i = 0; i < parNumber; i++) {

        auto &sphere = systemSP[i];
        bool swimmerHead = !isOdd(sphere.gid);

        auto &springCenter = dataFinder.findData[i];
        const auto Qv = springCenter - sphere.pos;
        const double Qabs = normVec3(Qv);
        PS::F64vec3 force;
        force = (tanh((1 / (ratio - 1)) * ((Qabs / l) / ratio - 1)) + 1) * Fmax * 0.5 * Qv / Qabs;
//        force = 0;    // (tanh((1 / (ratio - 1)) * ((Qabs / l) / ratio - 1)) + 1) * Fmax * 0.5 * Qv / Qabs;
//        force = ((Qabs - l) * h) * Qv;

        if (Qabs >= l * 1.5) {
            // rare event, larger than the max length, the FENE formula gives repulsive force
            printf("------------------------------------------------\n");
            printf("spring length larger than 1.5 L0 : %f\n", Qabs);
            printf("%f,%f,%f\n", sphere.pos.x, sphere.pos.y, sphere.pos.z);
            printf("%f,%f,%f\n", sphere.forceCollision.x, sphere.forceCollision.y, sphere.forceCollision.z);
            printf("%f,%f,%f\n", springCenter.x, springCenter.y, springCenter.z);
            printf("------------------------------------------------\n");
            lmN += 1;
        } else {
//             force = (h * (1 - l / Qabs) / (1 - pow((Qabs - l) / (lm - l), 2))) * Qv;
        }
        sphere.forceNear = force;
        if (swimmerHead) { // swimmer force
            sphere.forceExt = -(runConfig.swimForce / Qabs) * Qv;
        } else {
            sphere.forceExt = 0;
        }

        zeroVec3(sphere.torqueNear);
        zeroVec3(sphere.torqueExt);
    }

    printf(" number of spring length larger than 1.5 L0 : %d\n", lmN / 2);

}

void SphereSystem::calcExtForce() {
}

//#pragma GCC pop_options
