#include <iostream>

#include "Sphere/Sphere.hpp"
#include "Util/Buffer.hpp"

void testSimple() {
    std::string name("SPHERE");
    std::array<double, 4> array4{1, 2, 3, 4};
    std::vector<char> buf;

    Buffer mybuffer(buf);
    mybuffer.pack(name);
    mybuffer.pack(array4);
    name.clear();
    array4.fill(0);

    mybuffer.setReadPos(0);
    mybuffer.unpack(name, buf);
    mybuffer.unpack(array4, buf);
    std::cout << name << " " << array4[0] << array4[1] << array4[2] << array4[3] << std::endl;
}

void testSingleSphere() {
    std::vector<char> buf;
    {
        Sphere sphere(5, 1.0, 1.5, Evec3(1, 2, 3), Equatn(1, 0, 0, 0));
        sphere.addLayer(std::string("testLayer"), Shexp::KIND::STK, 4);
        auto &specValue = sphere.sphLayer["testLayer"]->gridValue;
        std::fill(specValue.begin(), specValue.end(), 1.0);
        sphere.Pack(buf);
    }
    Sphere sphere2;
    sphere2.Unpack(buf);
    sphere2.dumpSphere();
    sphere2.dumpLayer(std::string("testLayer"));
}

void testMultiSphere() {
    std::vector<char> buf1;
    std::vector<char> buf2;
    { // pack sphere 1
        Sphere sphere(5, 1.0, 1.5, Evec3(1, 2, 3), Equatn(1, 0, 0, 0));
        sphere.addLayer(std::string("testLayer"), Shexp::KIND::STK, 4);
        auto &specValue = sphere.sphLayer["testLayer"]->gridValue;
        std::fill(specValue.begin(), specValue.end(), 1.0);
        sphere.Pack(buf1);
    }
    { // pack sphere 2
        Sphere sphere(8, 2.0, 3.5, Evec3(3, 2, 1), Equatn(0, 1, 0, 0));
        sphere.addLayer(std::string("testLayer"), Shexp::KIND::LAP, 3);
        auto &specValue = sphere.sphLayer["testLayer"]->gridValue;
        std::fill(specValue.begin(), specValue.end(), 2.0);
        sphere.Pack(buf2);
    }
    Sphere sphere1;
    sphere1.Unpack(buf1);
    sphere1.dumpSphere();
    sphere1.dumpLayer(std::string("testLayer"));
    Sphere sphere2;
    sphere2.Unpack(buf2);
    sphere2.dumpSphere();
    sphere2.dumpLayer(std::string("testLayer"));
}

int main() {
    testSimple();
    testSingleSphere();
    testMultiSphere();
    return 0;
}