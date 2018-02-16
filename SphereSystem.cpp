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
    id_snap = 0;

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

}