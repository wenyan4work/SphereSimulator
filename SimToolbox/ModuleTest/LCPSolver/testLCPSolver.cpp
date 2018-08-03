#include "../../Collision/CPSolver.hpp"
#include <mpi.h>

#include <vector>

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    const int globalSize = argc > 1 ? atoi(argv[1]) : 500;
    const double diagonal = argc > 2 ? atof(argv[2]) : 1.0;
    CPSolver lcptest(globalSize, diagonal);

    double tol = 1e-5;

    // // mmNewton
    // lcptest.test_LCP(tol, globalSize, 0);
    // // APGD
    // lcptest.test_LCP(tol, globalSize, 1);
    // BBPGD
    lcptest.test_LCP(tol, globalSize, 2);
    // // APGD+Newton
    // lcptest.test_LCP(tol, globalSize, 3);
    // // BBPGD+mmNewton
    // lcptest.test_LCP(tol, globalSize, 4);

    MPI_Finalize();
    return 0;
}