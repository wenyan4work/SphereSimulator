
#include <cstdio>
#include <vector>

#include "Sylinder.hpp"
#include "Util/Base64.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/IOHelper.hpp"

/*****************************************************
 *  Sphero-cylinder
 ******************************************************/

Sylinder::Sylinder(const int &gid_, const double &radius_, const double &radiusCollision_, const double &length_,
                   const double &lengthCollision_, const Evec3 &pos_, const Equatn &orientation_) {
    gid = gid_;
    radius = radius_;
    radiusCollision = radiusCollision_;
    length = length_;
    lengthCollision = lengthCollision_;
    pos = pos_;
    orientation = orientation_;

    clear();
    return;
}

void Sylinder::clear() {
    vel.setZero();
    omega.setZero();
    velCol.setZero();
    omegaCol.setZero();
    velBrown.setZero();
    omegaBrown.setZero();
    velHydro.setZero();
    omegaHydro.setZero();
}

void Sylinder::dumpSylinder() const {
    printf("gid %8d, r %8f, rCol %8f, l %8f, lCol %8f, pos %8f, %8f, %8f\n", gid, radius, radiusCollision, length,
           lengthCollision, pos[0], pos[1], pos[2]);
    printf("vel %8f, %8f, %8f; omega %8f, %8f, %8f\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %8f, %8f, %8f, %8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
}

void Sylinder::Pack(std::vector<char> &buff) const {
    Buffer buffer(buff);
    // head
    buffer.pack(std::string("SYLINDER"));
    // fixed size data
    buffer.pack(gid);                                                 // int gid = INVALID;
    buffer.pack(globalIndex);                                         // int gid = INVALID;
    buffer.pack(radius);                                              // double radius;
    buffer.pack(radiusCollision);                                     // double radiusCollision;
    buffer.pack(length);                                              // double ;
    buffer.pack(lengthCollision);                                     // double ;
    buffer.pack(radiusSearch);                                        // double
    buffer.pack(std::array<double, 3>{pos[0], pos[1], pos[2]});       // Evec3 pos;
    buffer.pack(std::array<double, 3>{vel[0], vel[1], vel[2]});       // Evec3 vel;
    buffer.pack(std::array<double, 3>{omega[0], omega[1], omega[2]}); // Evec3 omega;
    buffer.pack(std::array<double, 4>{orientation.w(), orientation.x(), orientation.y(),
                                      orientation.z()}); // Equatn orientation;
}

void Sylinder::Unpack(const std::vector<char> &buff) {
    Buffer buffer;
    // head
    std::string strbuf;
    buffer.unpack(strbuf, buff);
    assert(strbuf == std::string("SYLINDER"));
    // fixed size data
    buffer.unpack(gid, buff);             // int gid = INVALID;
    buffer.unpack(globalIndex, buff);     // int gid = INVALID;
    buffer.unpack(radius, buff);          // double radius;
    buffer.unpack(radiusCollision, buff); // double radiusCollision;
    buffer.unpack(length, buff);          // double
    buffer.unpack(lengthCollision, buff); // double
    buffer.unpack(radiusSearch, buff);    // double
    std::array<double, 3> array3;
    buffer.unpack(array3, buff); // Evec3 pos;
    pos[0] = array3[0];
    pos[1] = array3[1];
    pos[2] = array3[2];
    buffer.unpack(array3, buff); // Evec3 vel;
    vel[0] = array3[0];
    vel[1] = array3[1];
    vel[2] = array3[2];
    buffer.unpack(array3, buff); // Evec3 omega;
    omega[0] = array3[0];
    omega[1] = array3[1];
    omega[2] = array3[2];
    std::array<double, 4> array4;
    buffer.unpack(array4, buff); // Equatn orientation;
    orientation.w() = array4[0];
    orientation.x() = array4[1];
    orientation.y() = array4[2];
    orientation.z() = array4[3];
}

void Sylinder::writeVTP(const std::vector<Sylinder> &sylinder, const std::string &prefix, const std::string &postfix,
                        int rank) {
    // for each sylinder:
    /*
     Procedure for dumping sylinders in the system:
    each sylinder writes a polyline with two (connected) points. points are labeled with float -1 and 1
    sylinder data written as cell data
    1 basic data (x,y,z), r, rcol, v, omega, etc output as VTK POLY_DATA file (*.vtp), each partiel is a VTK_VERTEX
    For each dataset, rank 0 writes the parallel header , then each rank write its own serial vtp/vtu file
    */

    // write VTP for basic data
    //  use float to save some space
    const int sylinderNumber = sylinder.size();
    // point and point data
    std::vector<double> pos(6 * sylinderNumber); // position always in Float64
    std::vector<float> label(2 * sylinderNumber);

    // point connectivity of line
    std::vector<int32_t> connectivity(2 * sylinderNumber);
    std::vector<int32_t> offset(sylinderNumber);

    // sylinder data
    std::vector<int> gid(sylinderNumber);
    std::vector<float> radius(sylinderNumber);
    std::vector<float> radiusCollision(sylinderNumber);
    std::vector<float> length(sylinderNumber);
    std::vector<float> lengthCollision(sylinderNumber);
    std::vector<float> vel(3 * sylinderNumber);
    std::vector<float> omega(3 * sylinderNumber);
    std::vector<float> velCol(3 * sylinderNumber);
    std::vector<float> omegaCol(3 * sylinderNumber);
    std::vector<float> velHydro(3 * sylinderNumber);
    std::vector<float> omegaHydro(3 * sylinderNumber);
    std::vector<float> velBrown(3 * sylinderNumber);
    std::vector<float> omegaBrown(3 * sylinderNumber);
    std::vector<float> xnorm(3 * sylinderNumber);
    std::vector<float> znorm(3 * sylinderNumber);

#pragma omp parallel for
    for (int i = 0; i < sylinderNumber; i++) {
        const auto &sy = sylinder[i];
        // point and point data
        Evec3 direction = sy.orientation * Evec3(0, 0, 1);
        Evec3 end0 = sy.pos - direction * (sy.length * 0.5);
        Evec3 end1 = sy.pos + direction * (sy.length * 0.5);
        pos[6 * i] = end0[0];
        pos[6 * i + 1] = end0[1];
        pos[6 * i + 2] = end0[2];
        pos[6 * i + 3] = end1[0];
        pos[6 * i + 4] = end1[1];
        pos[6 * i + 5] = end1[2];
        label[2 * i] = -1;
        label[2 * i + 1] = 1;

        // connectivity
        connectivity[2 * i] = 2 * i;         // index of point 0 in line
        connectivity[2 * i + 1] = 2 * i + 1; // index of point 1 in line
        offset[i] = 2 * i + 2;               // offset is the end of each line. in fortran indexing

        // sylinder data
        gid[i] = sy.gid;
        radius[i] = sy.radius;
        radiusCollision[i] = sy.radiusCollision;
        length[i] = sy.length;
        lengthCollision[i] = sy.lengthCollision;

        Evec3 nx = sy.orientation * Evec3(1, 0, 0);
        Evec3 nz = sy.orientation * Evec3(0, 0, 1);
        for (int j = 0; j < 3; j++) {
            vel[3 * i + j] = sy.vel[j];
            omega[3 * i + j] = sy.omega[j];
            velBrown[3 * i + j] = sy.velBrown[j];
            omegaBrown[3 * i + j] = sy.omegaBrown[j];
            velCol[3 * i + j] = sy.velCol[j];
            omegaCol[3 * i + j] = sy.omegaCol[j];
            velHydro[3 * i + j] = sy.velHydro[j];
            omegaHydro[3 * i + j] = sy.omegaHydro[j];

            xnorm[3 * i + j] = nx[j];
            znorm[3 * i + j] = nz[j];
        }
    }

    std::ofstream file(prefix + std::string("Sylinder_") + "r" + std::to_string(rank) + std::string("_") + postfix +
                           std::string(".vtp"),
                       std::ios::out);

    IOHelper::writeHeadVTP(file);

    file << "<Piece NumberOfPoints=\"" << sylinderNumber * 2 << "\" NumberOfLines=\"" << sylinderNumber << "\">\n";
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
    file << "</PointData>\n";
    // cell data
    file << "<CellData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
    IOHelper::writeDataArrayBase64(radius, "radius", 1, file);
    IOHelper::writeDataArrayBase64(radiusCollision, "radiusCollision", 1, file);
    IOHelper::writeDataArrayBase64(length, "length", 1, file);
    IOHelper::writeDataArrayBase64(lengthCollision, "lengthCollision", 1, file);
    IOHelper::writeDataArrayBase64(vel, "velocity", 3, file);
    IOHelper::writeDataArrayBase64(omega, "omega", 3, file);
    IOHelper::writeDataArrayBase64(velBrown, "velocityBrown", 3, file);
    IOHelper::writeDataArrayBase64(omegaBrown, "omegaBrown", 3, file);
    IOHelper::writeDataArrayBase64(velCol, "velocityCollision", 3, file);
    IOHelper::writeDataArrayBase64(omegaCol, "omegaCollision", 3, file);
    IOHelper::writeDataArrayBase64(velHydro, "velocityHydro", 3, file);
    IOHelper::writeDataArrayBase64(omegaHydro, "omegaHydro", 3, file);
    IOHelper::writeDataArrayBase64(xnorm, "xnorm", 3, file);
    IOHelper::writeDataArrayBase64(znorm, "znorm", 3, file);
    file << "</CellData>\n";
    file << "</Piece>\n";

    IOHelper::writeTailVTP(file);
    file.close();
}

void Sylinder::writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs) {
    std::vector<std::string> pieceNames;

    std::vector<IOHelper::FieldVTU> pointDataFields;
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "endLabel");

    std::vector<IOHelper::FieldVTU> cellDataFields;
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radius");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radiusCollision");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "length");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "lengthCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocity");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omega");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityBrown");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaBrown");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaCollision");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "velocityHydro");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omegaHydro");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "xnorm");
    cellDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "znorm");

    for (int i = 0; i < nProcs; i++) {
        pieceNames.emplace_back(std::string("Sylinder_") + std::string("r") + std::to_string(i) + "_" + postfix +
                                ".vtp");
    }

    IOHelper::writePVTPFile(prefix + "Sylinder_" + postfix + ".pvtp", pointDataFields, cellDataFields, pieceNames);
}

void Sylinder::stepEuler(double dt) {
    pos += vel * dt;
    EquatnHelper::rotateEquatn(orientation, omega, dt);
}