#include "../../SPHExp.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

int main() {

    SPHExp sph(SPHExp::KIND::LAP, "test", 12);
    sph.spectralCoeff[0] = 1;
    sph.spectralCoeff[1] = 2;

    std::ofstream vtkfile("test.vtu", std::ios::out);

    vtkfile << "<?xml version=\"1.0\"?>\n";
    vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"  header_type=\"UInt32\">\n";
    vtkfile << "<UnstructuredGrid>\n";
    sph.dumpVTK(vtkfile);
    vtkfile << "</UnstructuredGrid>\n";
    vtkfile << "</VTKFile>" << std::endl;

    vtkfile.close();

    return 0;
}
