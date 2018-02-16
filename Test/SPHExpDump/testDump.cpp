
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "Sphere/SPHExp.hpp"
#include "Sphere/Sphere.hpp"
#include "Util/IOHelper.hpp"

void testVTK() {
    std::vector<Sphere> sphere(10);
    for (int i = 0; i < 10; i++) {
        // fill random data
        auto &s = sphere[i];
        s.gid = i;
        s.globalIndex = i;
        s.radius = drand48();
        s.radiusCollision = drand48();
        Emap3(s.pos).setRandom();
        Emap3(s.vel).setRandom();
        Emap3(s.omega).setRandom();
        s.orientation = Equatn::UnitRandom();
        s.addLayer("lap", SPHExp::KIND::LAP, 6, Equatn::UnitRandom());
        s.addLayer("stk", SPHExp::KIND::STK, 8, Equatn::UnitRandom());
    }
    Sphere::writeVTP(sphere, "000", 0);
    Sphere::writeVTU(sphere, "000", 0);
    Sphere::writePVTP("000", 1);
    std::vector<std::pair<int, std::string>> dataFields;
    std::vector<IOHelper::IOTYPE> types = {IOHelper::IOTYPE::Float64, IOHelper::IOTYPE::Float64,
                                           IOHelper::IOTYPE::Float64};
    dataFields.emplace_back(std::pair<int, std::string>(1, "lap"));
    dataFields.emplace_back(std::pair<int, std::string>(3, "stk"));
    Sphere::writePVTU(dataFields, types, "000", 1);
}

int main() {

    testVTK();
    return 0;
}

/*
 Procedure for dumping spheres in the system:

    1 basic data (x,y,z), r, rcol, v, omega, etc output as VTK POLY_DATA file (*.vtp), each partiel is a VTK_VERTEX
    2 each sph with the same name, output as UnstructuredGrid (*.vtu). each sph belongs to one particle is an
   independent piece.

    For each dataset, rank 0 writes the parallel header , then each rank write its own serial vtp/vtu file

*/