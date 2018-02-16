#include <cassert>
#include <cmath>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "SphereSystem.hpp"
#include "Util/IOHelper.hpp"

SphereSystem::SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv)
    : runConfig(configFile) {
    int mpiflag;
    MPI_Initialized(&mpiflag);
    assert(mpiflag);

    commRcp = getMPIWORLDTCOMM();

    sysTime = 0;
    snapTime = 0;
    snapID = 0;

    // TRNG pool must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(runConfig.rngSeed);

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

    // initialize force calculate tree
    // initial resolve and reset move to zero
    // 5 pass to resolve initial overlaps

    printf("SphereSystem initialized on process: %d \n", commRcp->getRank());
}

void SphereSystem::setInitial(const std::string &initPos) {}

void SphereSystem::writeVTK() {
    Sphere::writeVTP(sphere, std::to_string(snapID), commRcp->getRank());
    Sphere::writeVTU(sphere, std::to_string(snapID), commRcp->getRank());
    // write parallel head
    if (commRcp->getSize() > 1 && commRcp->getRank() == 0) {
        Sphere::writePVTP("000", commRcp->getSize());
        std::vector<std::pair<int, std::string>> dataFields;
        std::vector<IOHelper::IOTYPE> types = {IOHelper::IOTYPE::Float64, IOHelper::IOTYPE::Float64};
        for (auto &layer : sphere[0].sphLayer) {
            dataFields.emplace_back(std::pair<int, std::string>(layer.second->dimension, layer.second->name));
        }
        dataFields.emplace_back(std::pair<int, std::string>(3, "stk"));
        Sphere::writePVTU(dataFields, types, std::to_string(snapID), commRcp->getSize());
    }
}