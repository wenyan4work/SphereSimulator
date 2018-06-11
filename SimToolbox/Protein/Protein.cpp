
#include "Protein.hpp"

#include <cstdio>
#include <vector>

#include "Util/Base64.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/IOHelper.hpp"

void Protein::writeVTP(const std::vector<Protein> &protein, const std::string &prefix, const std::string &postfix,
                       int rank) {
    // each protein is a vtk polyline with two vertices and one cell

    // write VTP for basic data
    //  use float to save some space
    const int proteinNumber = protein.size();
    // point and point data
    std::vector<double> pos(6 * proteinNumber);  // position always in Float64
    std::vector<float> label(2 * proteinNumber); // label is for end 0 and end 1

    // point connectivity of line
    std::vector<int32_t> connectivity(2 * proteinNumber);
    std::vector<int32_t> offset(proteinNumber);

    // protein data
    std::vector<int32_t> gid(proteinNumber);
    std::vector<int32_t> idBind(2 * proteinNumber);

#pragma omp parallel for
    for (int i = 0; i < proteinNumber; i++) {
        const auto &p = protein[i];
        // point and point data
        const Evec3 &end0 = p.posEnd[0];
        const Evec3 &end1 = p.posEnd[1];
        pos[6 * i] = end0[0];
        pos[6 * i + 1] = end0[1];
        pos[6 * i + 2] = end0[2];
        pos[6 * i + 3] = end1[0];
        pos[6 * i + 4] = end1[1];
        pos[6 * i + 5] = end1[2];
        label[2 * i] = 0;
        label[2 * i + 1] = 1;

        // connectivity
        connectivity[2 * i] = 2 * i;         // index of point 0 in line
        connectivity[2 * i + 1] = 2 * i + 1; // index of point 1 in line
        offset[i] = 2 * i + 2;               // offset is the end of each line. in fortran indexing

        // protein data
        // point data
        idBind[2 * i] = p.idBind[0];
        idBind[2 * i + 1] = p.idBind[1];
        // cell data
        gid[i] = p.gid;
    }

    std::ofstream file(prefix + std::string("Protein_") + "r" + std::to_string(rank) + std::string("_") + postfix +
                           std::string(".vtp"),
                       std::ios::out);

    IOHelper::writeHeadVTP(file);
    std::string contentB64; // data converted to base64 format
    contentB64.reserve(2 * proteinNumber);

    file << "<Piece NumberOfPoints=\"" << proteinNumber * 2 << "\" NumberOfLines=\"" << proteinNumber << "\">\n";
    // Points
    file << "<Points>\n";
    IOHelper::writeDataArrayBase64(pos, "position", 3, file);
    file << "</Points>\n";
    // cell definition
    file << "<Lines>\n";
    IOHelper::writeDataArrayBase64(connectivity, "connectivity", 1, file);
    IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
    file << "</Lines>\n";
    // point data
    file << "<PointData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(label, "endLabel", 1, file);
    IOHelper::writeDataArrayBase64(idBind, "idBind", 1, file);
    file << "</PointData>\n";
    // cell data
    file << "<CellData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
    file << "</CellData>\n";
    file << "</Piece>\n";

    IOHelper::writeTailVTP(file);
    file.close();
}

void Protein::writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs) {
    std::vector<std::string> pieceNames;

    std::vector<IOHelper::FieldVTU> pointDataFields;
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "endLabel");
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "idBind");

    std::vector<IOHelper::FieldVTU> cellDataFields;
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");

    for (int i = 0; i < nProcs; i++) {
        pieceNames.emplace_back(std::string("Protein_") + std::string("r") + std::to_string(i) + "_" + postfix +
                                ".vtp");
    }

    IOHelper::writePVTPFile(prefix + "Protein_" + postfix + ".pvtp", pointDataFields, cellDataFields, pieceNames);
}