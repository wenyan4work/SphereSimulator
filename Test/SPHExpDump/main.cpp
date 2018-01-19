#include "../../SPHExp.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

void writeHeadVTU(std::ofstream &vtkfile) {
    vtkfile << "<?xml version=\"1.0\"?>\n";
    vtkfile
        << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"  header_type=\"UInt32\">\n";
    vtkfile << "<UnstructuredGrid>\n";
}

void writeTailVTU(std::ofstream &vtkfile) {
    vtkfile << "</UnstructuredGrid>\n";
    vtkfile << "</VTKFile>" << std::endl;
}

void testSPHExp2VTK() {
    SPHExp sph1(SPHExp::KIND::LAP, "test", 72);
    SPHExp sph2(SPHExp::KIND::LAP, "test", 36);

    std::ofstream vtkfile("test.vtu", std::ios::out);
    std::array<double, 3> coordBase = {2, 2, 2};

    writeHeadVTU(vtkfile);
    sph1.dumpVTK(vtkfile);
    sph2.dumpVTK(vtkfile, coordBase);
    writeTailVTU(vtkfile);

    vtkfile.close();
}

int main() {
    testSPHExp2VTK();
    return 0;
}
