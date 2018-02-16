#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <mpi.h>

#include "Util/EigenDef.hpp"
#include "SphereSystem.hpp"

int main(int argc, char **argv) {
    Eigen::initParallel();
    Eigen::setNbThreads(1); // disable threading in eigen

    MPI_Init(&argc, &argv);

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
            mySystem.stepEuler();
            t += mySystem.runConfig.dt;
            iStep++;
        }
        MPI_Barrier(MPI_COMM_WORLD); // barrier before destructing mySystem;
        printf("mySystem destroyed");
    }

    // mpi finalize
    // let the root rank wait for other
    MPI_Finalize();
    return 0;
}
