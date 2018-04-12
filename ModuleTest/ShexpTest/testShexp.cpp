#include "../../Sphere/Shexp.hpp"

#include "PointDistribution.h"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <random>

#include <mpi.h>

void testLAPConvert(const int order = 12, const int repeat = 1000) {
    printf("testing LAP Conversion %d times for order=%d\n", repeat, order);
    Shexp sh(Shexp::KIND::LAP, "test", order, 2.0, Equatn::UnitRandom());

    // set random grid values
    std::vector<double> gridValueInitial(sh.getGridDOF(), 0);
    randomUniformFill(gridValueInitial, -2, 2);

    sh.gridValue = gridValueInitial;

    std::vector<double> spectralCoeff(sh.getSpectralDOF(), 0);
    sh.calcSpectralCoeff(spectralCoeff.data());
    sh.calcGridValue(spectralCoeff.data(), gridValueInitial.data()); // remove the nullspace of random grid value

    // transform 1000 times
    for (int i = 0; i < repeat; i++) {
        // grid to spectral
        sh.calcSpectralCoeff(spectralCoeff.data());
        // spectral to grid
        sh.calcGridValue(spectralCoeff.data());
    }

    // check error
    checkError(gridValueInitial, sh.gridValue);

    for (int i = 0; i < gridValueInitial.size(); i++) {
        printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValue[i], gridValueInitial[i],
               sh.gridValue[i] - gridValueInitial[i]);
    }

    return;
}

void testSTKConvert(const int order = 12, const int repeat = 1000) {
    printf("testing STK Conversion %d times for order=%d\n", repeat, order);
    Shexp sh(Shexp::KIND::STK, "test", order, 2.0, Equatn::UnitRandom());

    // set random grid values
    std::vector<double> gridValueInitial(sh.getGridDOF(), 0);
    randomUniformFill(gridValueInitial, -2, 2);

    sh.gridValue = gridValueInitial;

    std::vector<double> spectralCoeff(sh.getSpectralDOF(), 0);
    sh.calcSpectralCoeff(spectralCoeff.data());
    sh.calcGridValue(spectralCoeff.data(), gridValueInitial.data()); // remove the nullspace of random grid value

    // transform 1000 times
    for (int i = 0; i < repeat; i++) {
        // grid to spectral
        sh.calcSpectralCoeff(spectralCoeff.data());
        // spectral to grid
        sh.calcGridValue(spectralCoeff.data());
    }

    // check error
    checkError(gridValueInitial, sh.gridValue);

    for (int i = 0; i < gridValueInitial.size(); i++) {
        printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValue[i], gridValueInitial[i],
               sh.gridValue[i] - gridValueInitial[i]);
    }

    return;
}

void testSTKSL() {}

void testSTKDL() {}

void testTrac() {}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    testLAPConvert(24, 1000);
    testSTKConvert(24, 1000);

    MPI_Finalize();
    return 0;
}