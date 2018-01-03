
#include "../../Buffer.hpp"
#include "../../Sphere.hpp"
#include <iostream>

void testSimple() {
    std::string name("SPHERE");
    std::array<double, 4> array4{1, 2, 3, 4};
    Buffer mybuffer;
    mybuffer.pack(name);
    mybuffer.pack(array4);
    name.clear();
    array4.fill(0);

    mybuffer.unpack(name);
    mybuffer.unpack(array4);
    std::cout << name << " " << array4[0] << array4[1] << array4[2] << array4[3] << std::endl;
}

void testSingleSphere() {
    Buffer mybuffer;
    {
        Sphere sphere(5, 1.0, 1.5, Evec3(1, 2, 3), Equatn(1, 0, 0, 0));
        sphere.addLayer(std::string("testLayer"), SPHExp::KIND::STKDL, 4);
        auto &specValue = sphere.sphLayer["testLayer"]->spectralCoeff;
        std::fill(specValue.begin(), specValue.end(), 1.0);
        sphere.pack(mybuffer);
    }
    Sphere sphere2;
    sphere2.unpack(mybuffer);
    sphere2.dumpSphere();
    sphere2.dumpLayer(std::string("testLayer"));
}

void testMultiSphere() {
    Buffer mybuffer;
    { // pack sphere 1
        Sphere sphere(5, 1.0, 1.5, Evec3(1, 2, 3), Equatn(1, 0, 0, 0));
        sphere.addLayer(std::string("testLayer"), SPHExp::KIND::STKDL, 4);
        auto &specValue = sphere.sphLayer["testLayer"]->spectralCoeff;
        std::fill(specValue.begin(), specValue.end(), 1.0);
        sphere.pack(mybuffer);
    }
    { // pack sphere 2
        Sphere sphere(8, 2.0, 3.5, Evec3(3, 2, 1), Equatn(0, 1, 0, 0));
        sphere.addLayer(std::string("testLayer"), SPHExp::KIND::LAPDL, 3);
        auto &specValue = sphere.sphLayer["testLayer"]->spectralCoeff;
        std::fill(specValue.begin(), specValue.end(), 2.0);
        sphere.pack(mybuffer);
    }
    Sphere sphere1;
    sphere1.unpack(mybuffer);
    sphere1.dumpSphere();
    sphere1.dumpLayer(std::string("testLayer"));
    Sphere sphere2;
    sphere2.unpack(mybuffer);
    sphere2.dumpSphere();
    sphere2.dumpLayer(std::string("testLayer"));
}

int main() {
    testSimple();
    testSingleSphere();
    testMultiSphere();
    return 0;
}