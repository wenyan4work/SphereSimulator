#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "EigenDef.hpp"
#include "mpi.h"
#include "SphereSystem.hpp"

int main(int argc, char **argv) {
    Eigen::initParallel();
    Eigen::setNbThreads(1); // disable threading in eigen

    MPI_Init(&argc, &argv);
    // int myRank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // create a system and distribute it to all ranks
    std::string configFile = "runConfig.txt";
    std::string posFile = "posInitial.txt";
    {
        SphereSystem mySystem(configFile, posFile, argc, argv); // MPI is initialized inside PS::Initialize()
        // main time loop
        double t = 0;
        int iStep = 0;

        while (t < mySystem.runConfig.timeTotal + mySystem.runConfig.dt / 2) {
            MPI_Barrier(MPI_COMM_WORLD); // barrier before destructing mySystem;
            mySystem.stepForward();
            t += mySystem.runConfig.dt;
            iStep++;
        }
        MPI_Barrier(MPI_COMM_WORLD); // barrier before destructing mySystem;
        printf("mySystem destroyed");
    }
    // clean up

    // mpi finalize
    // let the root rank wait for other
    MPI_Finalize();
    return 0;
}
