#include "../../Sphere/Shexp.hpp"

#include "PointDistribution.h"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <random>

#include <mpi.h>

void testLAPConvert(const int order = 12) {
    printf("testing LAP Conversion 1000 times for order=%d\n", order);
    Shexp sh(Shexp::KIND::LAP, "test", order, Equatn::Identity());

    // set random grid values
    std::vector<double> gridValueInitial(sh.getGridDOF(), 0);
    randomUniformFill(gridValueInitial, -2, 2);

    sh.gridValues = gridValueInitial;

    std::vector<double> spectralCoeff(sh.getSpectralDOF(), 0);
    sh.calcSpectralValue(spectralCoeff.data());
    sh.calcGridValue(spectralCoeff.data(), gridValueInitial.data()); // remove the nullspace of random grid value

    // transform 1000 times
    for (int i = 0; i < 1000; i++) {
        // grid to spectral
        sh.calcSpectralValue(spectralCoeff.data());
        // spectral to grid
        sh.calcGridValue(spectralCoeff.data());
    }

    // check error
    checkError(gridValueInitial, sh.gridValues);

    for (int i = 0; i < gridValueInitial.size(); i++) {
        printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValues[i], gridValueInitial[i],
               sh.gridValues[i] - gridValueInitial[i]);
    }

    return;
}

void testSTKConvert(const int order = 12) {
    printf("testing STK Conversion 1000 times for order=%d\n", order);
    Shexp sh(Shexp::KIND::STK, "test", 12, Equatn::Identity());

    // set random grid values
    std::vector<double> gridValueInitial(sh.getGridDOF(), 0);
    randomUniformFill(gridValueInitial, -2, 2);

    sh.gridValues = gridValueInitial;

    std::vector<double> spectralCoeff(sh.getSpectralDOF(), 0);
    sh.calcSpectralValue(spectralCoeff.data());
    sh.calcGridValue(spectralCoeff.data(), gridValueInitial.data()); // remove the nullspace of random grid value

    // transform 1000 times
    for (int i = 0; i < 1000; i++) {
        // grid to spectral
        sh.calcSpectralValue(spectralCoeff.data());
        // spectral to grid
        sh.calcGridValue(spectralCoeff.data());
    }

    // check error
    checkError(gridValueInitial, sh.gridValues);

    for (int i = 0; i < gridValueInitial.size(); i++) {
        printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValues[i], gridValueInitial[i],
               sh.gridValues[i] - gridValueInitial[i]);
    }

    return;
}

void testSTKSL() {}

void testSTKDL() {}

void testTrac() {}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    testLAPConvert();

    MPI_Finalize();
    return 0;
}