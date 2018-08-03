#include "../../Collision/CPSolver.hpp"
#include <mpi.h>

#include <vector>

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    const int globalSize = argc > 1 ? atoi(argv[1]) : 500;
    const double diagonal = argc > 2 ? atof(argv[2]) : 0.0;
    CPSolver lcptest(globalSize, diagonal);

    double tol = 1e-5;

    lcptest.test_LCP(tol, globalSize, 0); // mmNewton
    lcptest.test_LCP(tol, globalSize, 1); // APGD
    lcptest.test_LCP(tol, globalSize, 2); // BBPGD
    lcptest.test_LCP(tol, globalSize, 3); // APGD+Newton
    lcptest.test_LCP(tol, globalSize, 4); // BBPGD+mmNewton

    MPI_Finalize();
    return 0;
}