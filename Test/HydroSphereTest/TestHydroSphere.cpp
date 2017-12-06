/*
 * main.cpp
 *
 *  Created on: Oct 4, 2016
 *      Author: wyan
 */

#include "../../HydroSphere.h"
#include "mpi.h"

#include <climits>
#include <iomanip>
#include <iostream>
#include <unordered_map>

using namespace HydroSphere;

void testSphere(const int nLocal) {

    std::array<double, 3> boxLow = { 0, 0, 0 };
    std::array<double, 3> boxHigh = { 2, 2, 2 };
    HydroSphereSystem pointSphereSys(1.0, FMM_Wrapper::PAXIS::PXYZ, 10, boxLow, boxHigh);

    pointSphereSys.pointSphere.resize(nLocal);

    for (int i = 0; i < nLocal; i++) {
        for (int j = 0; j < 3; j++) {
            pointSphereSys.pointSphere[i].pos[j] = drand48() * (boxHigh[j] - boxLow[j]) + boxLow[j];
            pointSphereSys.pointSphere[i].radius = 1.0;
        }
        pointSphereSys.pointSphere[i].force[0] = 1.0;
    }

    pointSphereSys.updateSystem();
    pointSphereSys.farfieldMobilityApply();

    for (int i = 0; i < nLocal; i++) {
        const auto & sphere = pointSphereSys.pointSphere[i];
        printf("%lf,%lf,%lf\n", sphere.vel[0], sphere.vel[1], sphere.vel[2]);
    }

    pointSphereSys.calcFlowOnGrid(0.01, "test", 0, 0, 0);

}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    testSphere(1000);

    MPI_Finalize();
    return 0;
}
