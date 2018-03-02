#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <mpi.h>

#include "SphereSystem.hpp"
#include "Util/EigenDef.hpp"

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
            mySystem.step();
            t += mySystem.runConfig.dt;
            iStep++;
            if (iStep % 10 == 0) {
                mySystem.partition();
            }
        }
    }
    // mpi finalize
    // let the root rank wait for other
    MPI_Finalize();
    return 0;
}
