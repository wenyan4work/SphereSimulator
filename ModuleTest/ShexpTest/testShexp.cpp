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

    // for (int i = 0; i < gridValueInitial.size(); i++) {
    //     printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValue[i], gridValueInitial[i],
    //            sh.gridValue[i] - gridValueInitial[i]);
    // }

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

    // for (int i = 0; i < gridValueInitial.size(); i++) {
    //     printf("%18.16lf\t\t%18.16lf\t\t%g\n", sh.gridValue[i], gridValueInitial[i],
    //            sh.gridValue[i] - gridValueInitial[i]);
    // }

    return;
}

void testSTKSL(const int order = 12, const int npts = 1000, bool interior = false) {

    // generate a random radius
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.2, 5.0);
    double rad = dis(gen);
    double rtrg;

    if (interior) {
        printf("Interior testing STK SL on %d points for order=%d, from sphere radius %lf\n", npts, order, rad);
        rtrg = rad * 0.2; // a far point
    } else {
        printf("Exterior testing STK SL on %d points for order=%d, from sphere radius %lf\n", npts, order, rad);
        rtrg = rad * 4; // a far point
    }

    Shexp sh(Shexp::KIND::STK, "test", order, rad, Equatn::UnitRandom());
    sh.randomFill();
    std::vector<double> spectralCoeff(sh.getSpectralDOF());
    sh.calcSpectralCoeff(spectralCoeff.data());

    // generate trgs for exterior
    std::vector<double> trgXYZ(3 * npts);
    std::vector<double> trgValue(3 * npts, 0);

    // exterior points

    for (int i = 0; i < npts; i++) {
        Evec3 pos = Evec3::Random();
        pos.normalize();
        pos = pos * rtrg;
        trgXYZ[3 * i] = pos[0];
        trgXYZ[3 * i + 1] = pos[1];
        trgXYZ[3 * i + 2] = pos[2];
    }

    std::vector<double> trgXYZGrid = trgXYZ;

    // spectral SL eval
    sh.calcSDLNF(spectralCoeff.data(), npts, trgXYZ.data(), trgValue.data(), interior, true);

    // grid SL eval
    const double fac8pi = 1 / (8 * 3.14159265358979323846);
    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    sh.getGridWithPole(gridPoints, gridWeights);
    const int Ngrid = sh.getGridNumber();
    std::vector<double> trgValueGrid(3 * npts, 0);
    for (int i = 0; i < npts; i++) {
        double tvx = 0, tvy = 0, tvz = 0;
        const double tx = trgXYZGrid[3 * i];
        const double ty = trgXYZGrid[3 * i + 1];
        const double tz = trgXYZGrid[3 * i + 2];
        for (int j = 0; j < Ngrid; j++) {
            // stokes single layer kernel
            const double lx = gridPoints[3 * (j + 1)];
            const double ly = gridPoints[3 * (j + 1) + 1];
            const double lz = gridPoints[3 * (j + 1) + 2];
            const double fx = sh.gridValue[3 * j] * gridWeights[j + 1];
            const double fy = sh.gridValue[3 * j + 1] * gridWeights[j + 1];
            const double fz = sh.gridValue[3 * j + 2] * gridWeights[j + 1];
            const double rx = (tx - lx);
            const double ry = (ty - ly);
            const double rz = (tz - lz);
            const double rnorm2 = rx * rx + ry * ry + rz * rz;
            const double rinv = 1 / sqrt(rnorm2);
            const double rinv3 = rinv * rinv * rinv;
            const double commonFac = (rx * fx + ry * fy + rz * fz);
            tvx += fx * rinv + commonFac * rx * rinv3;
            tvy += fy * rinv + commonFac * ry * rinv3;
            tvz += fz * rinv + commonFac * rz * rinv3;
        }
        trgValueGrid[3 * i] = tvx * fac8pi;
        trgValueGrid[3 * i + 1] = tvy * fac8pi;
        trgValueGrid[3 * i + 2] = tvz * fac8pi;
    }

    // check error
    checkError(trgValue, trgValueGrid);

    // for (int i = 0; i < 3 * npts; i++) {
    //     printf("%18.16lf\t\t%18.16lf\t\t%g\n", trgValue[i], trgValueGrid[i], trgValue[i] - trgValueGrid[i]);
    // }

    return;
}

void testSTKDL(const int order = 12, const int npts = 1000, bool interior = false) {

    // generate a random radius
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.2, 5.0);
    double rad = dis(gen);
    double rtrg;

    if (interior) {
        printf("Interior testing STK DL on %d points for order=%d, from sphere radius %lf\n", npts, order, rad);
        rtrg = rad * 0.2; // a far point
    } else {
        printf("Exterior testing STK DL on %d points for order=%d, from sphere radius %lf\n", npts, order, rad);
        rtrg = rad * 4; // a far point
    }

    Shexp sh(Shexp::KIND::STK, "test", order, rad, Equatn::UnitRandom());
    sh.randomFill();
    std::vector<double> spectralCoeff(sh.getSpectralDOF());
    sh.calcSpectralCoeff(spectralCoeff.data());

    // generate trgs for exterior
    std::vector<double> trgXYZ(3 * npts);
    std::vector<double> trgValue(3 * npts, 0);

    // exterior points

    for (int i = 0; i < npts; i++) {
        Evec3 pos = Evec3::Random();
        pos.normalize();
        pos = pos * rtrg;
        trgXYZ[3 * i] = pos[0];
        trgXYZ[3 * i + 1] = pos[1];
        trgXYZ[3 * i + 2] = pos[2];
    }

    std::vector<double> trgXYZGrid = trgXYZ;

    // spectral DL eval
    sh.calcSDLNF(spectralCoeff.data(), npts, trgXYZ.data(), trgValue.data(), interior, false);

    // grid DL eval
    const double fac4pi = -3 / (4 * 3.14159265358979323846);
    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    std::vector<double> gridNorms;
    sh.getGridWithPole(gridPoints, gridWeights, Evec3::Zero(), &gridNorms);
    const int Ngrid = sh.getGridNumber();
    std::vector<double> trgValueGrid(3 * npts, 0);
    for (int i = 0; i < npts; i++) {
        double tvx = 0, tvy = 0, tvz = 0;
        const double tx = trgXYZGrid[3 * i];
        const double ty = trgXYZGrid[3 * i + 1];
        const double tz = trgXYZGrid[3 * i + 2];
        for (int j = 0; j < Ngrid; j++) {
            // stokes single layer kernel
            const double lx = gridPoints[3 * (j + 1)];
            const double ly = gridPoints[3 * (j + 1) + 1];
            const double lz = gridPoints[3 * (j + 1) + 2];
            const double fx = sh.gridValue[3 * j] * gridWeights[j + 1];
            const double fy = sh.gridValue[3 * j + 1] * gridWeights[j + 1];
            const double fz = sh.gridValue[3 * j + 2] * gridWeights[j + 1];
            const double nx = gridNorms[3 * (j + 1)]; // norm toward outside
            const double ny = gridNorms[3 * (j + 1) + 1];
            const double nz = gridNorms[3 * (j + 1) + 2];
            const double rx = (tx - lx);
            const double ry = (ty - ly);
            const double rz = (tz - lz);
            const double rnorm2 = rx * rx + ry * ry + rz * rz;
            const double rinv = 1 / sqrt(rnorm2);
            const double rinv3 = rinv * rinv * rinv;
            const double rinv5 = rinv3 * rinv * rinv;
            const double commonFac = rx * rx * fx * nx + rx * ry * fx * ny + rx * rz * fx * nz + ry * rx * fy * nx +
                                     ry * ry * fy * ny + ry * rz * fy * nz + rz * rx * fz * nx + rz * ry * fz * ny +
                                     rz * rz * fz * nz;
            tvx += rx * commonFac * rinv5;
            tvy += ry * commonFac * rinv5;
            tvz += rz * commonFac * rinv5;
        }
        trgValueGrid[3 * i] = tvx * fac4pi;
        trgValueGrid[3 * i + 1] = tvy * fac4pi;
        trgValueGrid[3 * i + 2] = tvz * fac4pi;
    }

    // check error
    checkError(trgValue, trgValueGrid);

    for (int i = 0; i < 3 * npts; i++) {
        printf("%18.16lf\t\t%18.16lf\t\t%g\n", trgValue[i], trgValueGrid[i], trgValue[i] - trgValueGrid[i]);
    }

    return;
}

void testTrac() {}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    // for Eigen::setRandom
    srand((unsigned int)time(0));

    const int order = 16;
    testSTKSL(order, 1000, false);
    testSTKSL(order, 1000, true);

    testSTKDL(order, 1000, false);
    testSTKDL(order, 1000, true);

    testLAPConvert(order, 1000);
    testSTKConvert(order, 1000);

    MPI_Finalize();
    return 0;
}