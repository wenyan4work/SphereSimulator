#ifndef IOHELPER_HPP
#define IOHELPER_HPP

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

// for filesystem
// #include <boost/filesystem.hpp>
#include <errno.h>
#include <sys/stat.h>

#include "Base64.hpp"

class IOHelper {
  public:
    enum class IOTYPE {
        UInt8,
        Int32,
        Float32,
        Float64

    };

    static std::string getTypeName(IOTYPE type) {
        std::string name;
        if (type == IOTYPE::UInt8) {
            name = "UInt8";
        } else if (type == IOTYPE::Int32) {
            name = "Int32";
        } else if (type == IOTYPE::Float32) {
            name = "Float32";
        } else if (type == IOTYPE::Float64) {
            name = "Float64";
        }
        return name;
    }

    static void makeSubFolder(const std::string folder) {

        // boost::filesystem::create_directories(folder);
        // if (!boost::filesystem::is_directory(folder)) {
        //     printf("Error creating directory!n");
        //     exit(1);
        // }

        const int err = mkdir(folder.c_str(), 0755); // make one folder at a time. parent folder must exist
        if (errno == EEXIST) {
            printf("Directory already exists.\n");
            return;
        }
        if (err != 0) {
            printf("errno: %d \n", errno);
            printf("Error creating directory!\n");
            exit(1);
        }
    }

    /*******************************
     * VTK unstructured grid  *
     ********************************/

    static void writeHeadVTU(std::ofstream &vtkfile) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<UnstructuredGrid>\n";
    }

    static void writeTailVTU(std::ofstream &vtkfile) {
        vtkfile << "</UnstructuredGrid>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }

    static void writePVTUFile(std::string filename, const std::vector<std::pair<int, std::string>> &dataFields,
                              const std::vector<IOTYPE> types, const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PUnstructuredGrid GhostLevel=\"0\"> \n";

        pvtufile << "<PPointData Scalars=\"scalars\">\n";

        for (int i = 0; i < dataFields.size(); i++) {
            auto &data = dataFields[i];
            // data.first = dimension
            // data.second = name
            std::string type = getTypeName(types[i]);
            pvtufile << "<PDataArray Name=\"" << data.second << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.first << "\" format=\"binary\"/>\n";
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

    /*******************************
     * VTK poly data  *
     ********************************/

    static void writeHeadVTP(std::ofstream &vtkfile) {
        vtkfile << "<?xml version=\"1.0\"?>\n";
        vtkfile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\"  "
                   "header_type=\"UInt32\">\n";
        vtkfile << "<PolyData>\n";
    }

    static void writeTailVTP(std::ofstream &vtkfile) {
        vtkfile << "</PolyData>\n";
        vtkfile << "</VTKFile>" << std::endl;
    }

    static void writePVTPFile(std::string filename, const std::vector<std::pair<int, std::string>> &dataFields,
                              const std::vector<IOTYPE> &types, const std::vector<std::string> &pieceNames) {

        std::ofstream pvtufile(filename, std::ios::out);

        pvtufile << "<?xml version=\"1.0\"?>\n";
        pvtufile << "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                    "header_type=\"UInt32\"> \n";
        pvtufile << "<PPolyData GhostLevel=\"0\"> \n";

        pvtufile << "<PPointData Scalars=\"scalars\">\n";

        for (int i = 0; i < dataFields.size(); i++) {
            auto &data = dataFields[i];
            // data.first = dimension
            // data.second = name
            std::string type = getTypeName(types[i]);
            pvtufile << "<PDataArray Name=\"" << data.second << "\" type=\"" << type << "\" NumberOfComponents=\""
                     << data.first << "\" format=\"binary\"/>\n";
        }

        pvtufile << "</PPointData>\n";

        pvtufile << "<PPoints> \n";
        pvtufile << "<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\"/>\n";
        pvtufile << "</PPoints> \n";

        for (const auto &piece : pieceNames) {

            pvtufile << "<Piece Source=\"" << piece << "\"/>\n";
        }
        pvtufile << "</PPolyData>\n";
        pvtufile << "</VTKFile>\n";
        pvtufile.close();
    }

    /*
     * general VTK helper
     */
    template <class T>
    static void writeDataArrayBase64(std::vector<T> &data, const std::string &name, int numComp, std::ofstream &file) {
        // set type name
        std::string vtktype;
        if (std::is_same<T, int>::value) {
            vtktype = "Int32";
        } else if (std::is_same<T, float>::value) {
            vtktype = "Float32";
        } else if (std::is_same<T, double>::value) {
            vtktype = "Float64";
        } else if (std::is_same<T, uint8_t>::value) {
            vtktype = "UInt8";
        }
        // radiusCollision
        std::string contentB64;
        contentB64.clear();
        B64Converter::getBase64FromVector(data, contentB64);

        file << "<DataArray Name=\"" << name << "\" type=\"" << vtktype << "\" NumberOfComponents=\"" << numComp
             << "\" format=\"binary\">\n";
        file << contentB64 << "\n";
        file << "</DataArray>\n";
    }
};

#endif