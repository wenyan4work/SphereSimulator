#include "Sphere.hpp"

#include "SimToolbox/Util/Base64.hpp"
#include "SimToolbox/Util/EquatnHelper.hpp"
#include "SimToolbox/Util/IOHelper.hpp"

#include <cstdio>
#include <vector>

/*****************************************************
 *  Sphere
 ******************************************************/

void swap(Sphere &sphA, Sphere &sphB) {
    using std::swap;
    swap(sphA.gid, sphB.gid);
    swap(sphA.globalIndex, sphB.globalIndex);
    swap(sphA.radius, sphB.radius);
    swap(sphA.radiusCollision, sphB.radiusCollision);
    copySwap(sphA.pos, sphB.pos);
    copySwap(sphA.vel, sphB.vel);
    copySwap(sphA.omega, sphB.omega);
    copySwap(sphA.orientation, sphB.orientation);

    // swap the containers
    sphA.sphLayer.swap(sphB.sphLayer);
    // sphA.sphNeighbor.swap(sphB.sphNeighbor);

    return;
}

// void SphereIO::dumpSphere() const {
//     printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", gid, radius, radiusCollision, pos[0], pos[1], pos[2]);
//     printf("vel %8f, %8f, %8f; omega %8f, %8f, %8f\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
//     printf("orient %8f, %8f, %8f, %8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
//     return;
// }

// Sphere::Sphere(const SphereIO &sphere)
//     : Sphere::Sphere(sphere.gid, sphere.radius, sphere.radiusCollision, sphere.pos, sphere.orientation) {
//     return;
// }

Sphere::Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const Evec3 &pos_,
               const Equatn &orientation_) noexcept
    : gid(gid_), radius(radius_), radiusCollision(radiusCollision_), pos(pos_), orientation(orientation_) {
    vel.setZero();
    omega.setZero();
    // sphNeighbor.reserve(10);
    return;
}

Sphere::~Sphere() noexcept {
    for (auto &l : sphLayer) {
        delete l.second;
        l.second = nullptr;
    }
    sphLayer.clear();
    // sphNeighbor.clear();
}

Sphere::Sphere(const Sphere &other) noexcept {
    gid = other.gid;
    radius = other.radius;
    radiusCollision = other.radiusCollision;
    pos = other.pos;
    vel = other.vel;
    omega = other.omega;
    orientation = other.orientation;

    // copy neighbors
    // sphNeighbor = other.sphNeighbor;

    // copy layers
    for (auto &l : other.sphLayer) {
        sphLayer[l.first] = new Shexp(*(l.second));
    }
}

Sphere::Sphere(Sphere &&other) noexcept {
    gid = other.gid;
    radius = other.radius;
    radiusCollision = other.radiusCollision;
    for (int i = 0; i < 3; i++) {
        pos[i] = other.pos[i];
        vel[i] = other.vel[i];
        omega[i] = other.omega[i];
    }
    orientation = other.orientation;

    // directly move for rvalue
    // sphNeighbor = std::move(other.sphNeighbor);
    sphLayer = std::move(other.sphLayer);
}

Sphere &Sphere::operator=(Sphere other) noexcept {
    using std::swap;
    swap(*this, other);
    return *this;
}

void Sphere::addLayer(const std::string &name_, const Shexp::KIND &kind_, const int &order_, const double &radius_,
                      const Equatn orientation_) {
    // the orientation of the layer can be different from the orientation of the sphere
    sphLayer[name_] = new Shexp(kind_, name_, order_, radius_ > 0 ? radius_ : radius, orientation_);
}

void Sphere::dumpSphere() const {
    printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", gid, radius, radiusCollision, pos[0], pos[1], pos[2]);
    printf("vel %8f, %8f, %8f; omega %8f, %8f, %8f\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %8f, %8f, %8f, %8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
}

// void Sphere::dumpNeighbor() const {
//     for (auto &nb : sphNeighbor) {
//         printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", nb.gid, nb.radius, nb.radiusCollision, nb.pos[0],
//                nb.pos[1], nb.pos[2]);
//         printf("Orientation: %8f, %8f, %8f, %8f; posRelative %8f, %8f, %8f\n", nb.orientation.w(),
//         nb.orientation.x(),
//                nb.orientation.y(), nb.orientation.z(), nb.posRelative[0], nb.posRelative[1], nb.posRelative[2]);
//     }
// }

void Sphere::dumpLayer(const std::string &name) const { sphLayer.find(name)->second->dumpGridValue(); }

void Sphere::Pack(std::vector<char> &buff) const {
    Buffer buffer(buff);
    // head
    buffer.pack(std::string("SPHERE"));
    // fixed size data
    buffer.pack(gid);                                                 // int gid = INVALID;
    buffer.pack(globalIndex);                                         // int gid = INVALID;
    buffer.pack(radius);                                              // double radius;
    buffer.pack(radiusCollision);                                     // double radiusCollision;
    buffer.pack(std::array<double, 3>{pos[0], pos[1], pos[2]});       // Evec3 pos;
    buffer.pack(std::array<double, 3>{vel[0], vel[1], vel[2]});       // Evec3 vel;
    buffer.pack(std::array<double, 3>{omega[0], omega[1], omega[2]}); // Evec3 omega;
    buffer.pack(std::array<double, 4>{orientation.w(), orientation.x(), orientation.y(),
                                      orientation.z()}); // Equatn orientation;
    // variable size data
    const int nLayer = sphLayer.size();
    buffer.pack(nLayer); // number of layers
    for (auto &layer : sphLayer) {
        // name
        buffer.pack(layer.first);
        // data
        buffer.pack(static_cast<int>(layer.second->kind)); //    const KIND kind;
        buffer.pack(layer.second->order);                  // const int order;        // the order of expansion p
        buffer.pack(layer.second->name);                   //    const std::string name; // the name of this quantity
        buffer.pack(layer.second->radius);                 //    const std::string name; // the name of this quantity
        const Equatn &orientation = layer.second->orientation;
        buffer.pack(std::array<double, 4>{orientation.w(), orientation.x(), orientation.y(),
                                          orientation.z()}); // Equatn orientation;
        buffer.pack(layer.second->gridValue);
    }
}

void Sphere::Unpack(const std::vector<char> &buff) {
    Buffer buffer;
    // head
    std::string strbuf;
    buffer.unpack(strbuf, buff);
    assert(strbuf == std::string("SPHERE"));
    // fixed size data
    buffer.unpack(gid, buff);             // int gid = INVALID;
    buffer.unpack(globalIndex, buff);     // int gid = INVALID;
    buffer.unpack(radius, buff);          // double radius;
    buffer.unpack(radiusCollision, buff); // double radiusCollision;
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

    // variable size data
    // sphNeighbor.clear();

    int nLayer = 0;
    buffer.unpack(nLayer, buff); // number of layers
    std::string strbuf2;
    sphLayer.clear();
    for (int i = 0; i < nLayer; i++) {
        buffer.unpack(strbuf, buff); // the name
        // data
        int kind, order;
        double radiusLayer;
        buffer.unpack(kind, buff);        //    const KIND kind;
        buffer.unpack(order, buff);       // const int order;        // the order of expansion p
        buffer.unpack(strbuf2, buff);     //    const std::string name; // the name of this quantity
        buffer.unpack(radiusLayer, buff); // radius
        buffer.unpack(array4, buff);      // Equatn orientation;
        Equatn orientTemp;
        orientTemp.w() = array4[0];
        orientTemp.x() = array4[1];
        orientTemp.y() = array4[2];
        orientTemp.z() = array4[3];
        // allocate
        sphLayer[strbuf] = new Shexp(static_cast<Shexp::KIND>(kind), strbuf2, order, radiusLayer, orientTemp);
        // unpack the sph values
        buffer.unpack(sphLayer[strbuf]->gridValue, buff);
    }
}

void Sphere::writeVTP(const std::vector<Sphere> &sphere, const std::string &prefix, const std::string &postfix,
                      int rank) {
    // for each sphere:
    /*
     Procedure for dumping spheres in the system:

    1 basic data (x,y,z), r, rcol, v, omega, etc output as VTK POLY_DATA file (*.vtp), each partiel is a VTK_VERTEX
    2 each sph with the same name, output as UnstructuredGrid (*.vtu). each sph belongs to one particle is an
   independent piece.
    For each dataset, rank 0 writes the parallel header , then each rank write its own serial vtp/vtu file
    */

    // write VTP for basic data
    //  use float to save some space
    const int sphereNumber = sphere.size();
    std::vector<int> gid(sphereNumber);
    std::vector<float> radius(sphereNumber);
    std::vector<float> radiusCollision(sphereNumber);
    std::vector<double> pos(3 * sphereNumber); // position always in Float64 format
    std::vector<float> vel(3 * sphereNumber);
    std::vector<float> omega(3 * sphereNumber);
    std::vector<float> xnorm(3 * sphereNumber);
    std::vector<float> znorm(3 * sphereNumber);
    std::vector<float> forceCol(3 * sphereNumber);
    std::vector<float> torqueCol(3 * sphereNumber);
    std::vector<float> forceNonCol(3 * sphereNumber);
    std::vector<float> torqueNonCol(3 * sphereNumber);

#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        gid[i] = sphere[i].gid;
        radius[i] = sphere[i].radius;
        radiusCollision[i] = sphere[i].radiusCollision;
        Evec3 nx = sphere[i].orientation * Evec3(1, 0, 0);
        Evec3 nz = sphere[i].orientation * Evec3(0, 0, 1);
        for (int j = 0; j < 3; j++) {
            pos[3 * i + j] = sphere[i].pos[j];
            vel[3 * i + j] = sphere[i].vel[j];
            omega[3 * i + j] = sphere[i].omega[j];
            xnorm[3 * i + j] = nx[j];
            znorm[3 * i + j] = nz[j];
            forceCol[3 * i + j] = sphere[i].forceCol[j];
            torqueCol[3 * i + j] = sphere[i].torqueCol[j];
            forceNonCol[3 * i + j] = sphere[i].forceNonCol[j];
            torqueNonCol[3 * i + j] = sphere[i].torqueNonCol[j];
        }
    }

    std::ofstream file(prefix + std::string("Sphere_") + "r" + std::to_string(rank) + std::string("_") + postfix +
                           std::string(".vtp"),
                       std::ios::out);

    IOHelper::writeHeadVTP(file);
    std::string contentB64; // data converted to base64 format
    contentB64.reserve(4 * sphereNumber);

    file << "<Piece NumberOfPoints=\"" << sphereNumber << "\" NumberOfCells=\"" << 0 << "\">\n";
    // point data
    file << "<PointData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
    IOHelper::writeDataArrayBase64(radius, "radius", 1, file);
    IOHelper::writeDataArrayBase64(radiusCollision, "radiusCollision", 1, file);
    IOHelper::writeDataArrayBase64(vel, "velocity", 3, file);
    IOHelper::writeDataArrayBase64(omega, "omega", 3, file);
    IOHelper::writeDataArrayBase64(xnorm, "xnorm", 3, file);
    IOHelper::writeDataArrayBase64(znorm, "znorm", 3, file);
    IOHelper::writeDataArrayBase64(forceCol, "forceCol", 3, file);
    IOHelper::writeDataArrayBase64(torqueCol, "torqueCol", 3, file);
    IOHelper::writeDataArrayBase64(forceNonCol, "forceNonCol", 3, file);
    IOHelper::writeDataArrayBase64(torqueNonCol, "torqueNonCol", 3, file);
    file << "</PointData>\n";
    // no cell data
    // Points
    file << "<Points>\n";
    IOHelper::writeDataArrayBase64(pos, "points", 3, file);
    file << "</Points>\n";
    file << "</Piece>\n";

    IOHelper::writeTailVTP(file);
    file.close();
}

void Sphere::writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs) {
    std::vector<IOHelper::FieldVTU> dataFields;
    dataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
    dataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radius");
    dataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "radiusCollision");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "vel");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "omega");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "xnorm");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "znorm");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "forceCol");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "torqueCol");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "forceNonCol");
    dataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "torqueNonCol");

    std::vector<std::string> pieceNames;
    for (int i = 0; i < nProcs; i++) {
        pieceNames.emplace_back(std::string("Sphere_") + std::string("r") + std::to_string(i) + "_" + postfix + ".vtp");
    }

    IOHelper::writePVTPFile(prefix + "Sphere_" + postfix + ".pvtp", dataFields, pieceNames);
}

void Sphere::writeVTU(const std::vector<Sphere> &sphere, const std::vector<IOHelper::FieldVTU> &dataFields,
                      const std::string &prefix, const std::string &postfix, int rank) {
    // each sphere has the same number and name of layers
    // dump one vtu file for each sph

    for (auto &data : dataFields) {
        // generate one independent vtu file for each dataField
        // data in Float64
        if (data.type != IOHelper::IOTYPE::Float64) {
            printf("VTU Type error, must be in float64\n");
            exit(1);
        }
        std::ofstream file(prefix + std::string("Sphere_") + data.name + std::string("_r") + std::to_string(rank) +
                               std::string("_") + postfix + std::string(".vtu"),
                           std::ios::out);
        IOHelper::writeHeadVTU(file);
        if (sphere.size() > 0) {
            for (auto &s : sphere) {
                const auto &it = s.sphLayer.find(data.name);
                Evec3 coordBase(s.pos[0], s.pos[1], s.pos[2]);
                if (it != s.sphLayer.end())
                    it->second->writeVTU(file, coordBase);
            }
        }
        IOHelper::writeTailVTU(file);
    }
}

void Sphere::writePVTU(const std::vector<IOHelper::FieldVTU> &dataFields, const std::string &prefix,
                       const std::string &postfix, const int nProcs) {

    for (auto &data : dataFields) {
        // also write weight for each data
        std::vector<IOHelper::FieldVTU> fields;
        fields.push_back(data);
        fields.emplace_back(1, IOHelper::IOTYPE::Float64, std::string("weight"));
        std::vector<std::string> pieceNames;
        for (int i = 0; i < nProcs; i++) {
            pieceNames.emplace_back("Sphere_" + data.name + std::string("_r") + std::to_string(i) + "_" + postfix +
                                    ".vtu");
        }
        IOHelper::writePVTUFile(prefix + "Sphere_" + data.name + "_" + postfix + ".pvtu", fields, pieceNames);
    }
}

void Sphere::stepEuler(double dt) {
    pos += vel * dt;
    EquatnHelper::rotateEquatn(orientation, omega, dt);
    for (auto &layer : this->sphLayer) {
        EquatnHelper::rotateEquatn(layer.second->orientation, omega, dt);
    }
}

Shexp &Sphere::getLayer(const std::string &name) { return *(sphLayer.find(name)->second); }

const Shexp &Sphere::getLayer(const std::string &name) const { return *(sphLayer.find(name)->second); }