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

void writePVTUFile(std::string filename, const std::vector<std::pair<int, std::string>> &dataFields,
                   const std::vector<std::string> &pieceNames) {

    std::ofstream pvtufile(filename, std::ios::out);

    pvtufile << "<?xml version=\"1.0\"?>\n";
    pvtufile
        << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\"> \n";
    pvtufile << "<PUnstructuredGrid GhostLevel=\"0\"> \n";

    pvtufile << "<PPointData Scalars=\"scalars\">\n";

    for (const auto &data : dataFields) {
        // data.first = dimension
        // data.second = name
        pvtufile << "<PDataArray Name=\"" << data.second
                 << "\" type=\"Float64\" NumberOfComponents=\"1\" format=\"binary\"/>\n";
    }

    pvtufile << "</PPointData>\n";

    pvtufile << "<PPoints> \n";
    pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
    pvtufile << "</PPoints> \n";

    for (const auto &piece : pieceNames) {

        pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
    }
    pvtufile << "</PUnstructuredGrid>\n";
    pvtufile << "</VTKFile>\n";
    pvtufile.close();
}

void testSPHExp2VTK() {
    SPHExp sph1(SPHExp::KIND::LAP, "test", 72);
    SPHExp sph2(SPHExp::KIND::LAP, "test", 36);

    std::ofstream vtkfile("test1.vtu", std::ios::out);
    std::array<double, 3> coordBase = {2, 2, 2};

    writeHeadVTU(vtkfile);
    sph1.dumpVTK(vtkfile, {2, 0, 0});
    sph2.dumpVTK(vtkfile, {0, 2, 0});
    writeTailVTU(vtkfile);

    vtkfile.close();
}

int main() {
    std::vector<std::pair<int, std::string>> dataFields;
    std::vector<std::string> pieceNames;
    dataFields.emplace_back(std::pair<int, std::string>(1, "test"));
    dataFields.emplace_back(std::pair<int, std::string>(3, "weight"));

    pieceNames.emplace_back("test0.vtu");
    pieceNames.emplace_back("test1.vtu");

    writePVTUFile("test.pvtu", dataFields, pieceNames);

    // testSPHExp2VTK();
    return 0;
}

/*
    TODO: Procedure for dumping spheres in the system:

    1 basic data (x,y,z), r, rcol, v, omega, etc output as VTK POLY_DATA file (*.vtp), each partiel is a VTK_VERTEX
    2 each sph with the same name, output as UnstructuredGrid (*.vtu). each sph belongs to one particle is an
   independent piece.

    For each dataset, rank 0 writes the parallel header , then each rank write its own serial vtp/vtu file

*/