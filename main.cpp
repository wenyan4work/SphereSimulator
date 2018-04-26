#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <mpi.h>

#include "SphereSystem.hpp"
#include "Util/EigenDef.hpp"

int main(int argc, char **argv)
{
    Eigen::initParallel();
    Eigen::setNbThreads(1); // disable threading in eigen

    double t_start, t_stop, duration, max_duration;
    int rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();

    // create a system and distribute it to all ranks
    std::string configFile = "runConfig.txt";
    std::string posFile = "posInitial.txt";
    {
        SphereSystem mySystem(configFile, posFile, argc, argv); // MPI is initialized inside PS::Initialize()
        // main time loop
        double t = 0;
        int iStep = 0;
        while (t < mySystem.runConfig.timeTotal + mySystem.runConfig.dt / 2)
        {
            mySystem.step();
            t += mySystem.runConfig.dt;
            iStep++;
            if (iStep % 10 == 0)
            {
                mySystem.partition();
            }
            if (iStep % mySystem.runConfig.snapFreq == 0 && rank == 0)
            {
                std::cout << "t = " << t;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t_stop = MPI_Wtime();

    duration = t_stop - t_start;
    MPI_Reduce(&duration, &max_duration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "Runtime: " << t_stop - t_start << " seconds" << std::endl;
    }

    MPI_Finalize();

    return 0;
}
