
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "Sphere/Shexp.hpp"
#include "Sphere/Sphere.hpp"
#include "Util/IOHelper.hpp"

void testVTK() {
    const int num = 200;
    int k = 200 / num;
    int low = k * num, high = k * num + num - 1;
    IOHelper::makeSubFolder("./result");
    std::string baseFolder =
        "./result/" + std::to_string(low) + std::string("-") + std::to_string(high) + std::string("/");
    IOHelper::makeSubFolder(baseFolder);

    std::vector<Sphere> sphere(3);
    for (int i = 0; i < sphere.size(); i++) {
        // fill random data
        auto &s = sphere[i];
        s.gid = i;
        s.globalIndex = i;
        s.radius = 2.0 * (drand48() + 0.4);
        s.radiusCollision = s.radius * 1.2;
        s.pos = Evec3::Random() * 5;
        s.vel.setRandom();
        s.omega.setRandom();
        s.orientation = Equatn::UnitRandom();
        s.addLayer("lap", Shexp::KIND::LAP, 16, -1, Equatn::UnitRandom());
        s.getLayer("lap").randomFill(i);
        s.addLayer("stk", Shexp::KIND::STK, 22, -1, Equatn::UnitRandom());
        s.getLayer("stk").randomFill(i);
    }
    // VTP files
    Sphere::writeVTP(sphere, baseFolder, "0", 0);
    Sphere::writePVTP(baseFolder, "0", 1);

    // VTU files
    std::vector<IOHelper::FieldVTU> dataFields;
    dataFields.emplace_back(1, IOHelper::IOTYPE::Float64, std::string("lap"));
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float64, std::string("stk"));
    Sphere::writeVTU(sphere, dataFields, baseFolder, "0", 0);
    Sphere::writePVTU(dataFields, baseFolder, "0", 1);
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