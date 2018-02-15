#include "Sphere.hpp"
#include <cstdio>

using Emap3 = Eigen::Map<Evec3>;

void swap(Sphere &sphA, Sphere &sphB) {
    using std::swap;
    swap(sphA.gid, sphB.gid);
    swap(sphA.radius, sphB.radius);
    swap(sphA.radiusCollision, sphB.radiusCollision);

    Eigen::Map<Evec3>(sphA.pos).swap(Eigen::Map<Evec3>(sphB.pos));
    Eigen::Map<Evec3>(sphA.vel).swap(Eigen::Map<Evec3>(sphB.vel));
    Eigen::Map<Evec3>(sphA.omega).swap(Eigen::Map<Evec3>(sphB.omega));

    Equatn otemp = sphB.orientation;
    sphB.orientation = sphA.orientation;
    sphA.orientation = otemp;

    // swap the containers
    sphA.sphLayer.swap(sphB.sphLayer);
    sphA.sphNeighbor.swap(sphB.sphNeighbor);

    return;
}

void SphereIO::dumpSphere() const {
    printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", gid, radius, radiusCollision, pos[0], pos[1], pos[2]);
    printf("vel %8f, %8f, %8f; omega %8f, %8f, %8f\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %8f, %8f, %8f, %8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
    return;
}

Sphere::Sphere(const SphereIO &sphere)
    : Sphere::Sphere(sphere.gid, sphere.radius, sphere.radiusCollision, sphere.pos, sphere.orientation) {
    return;
}

Sphere::Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const Evec3 &pos_,
               const Equatn &orientation_) noexcept
    : gid(gid_), radius(radius_), radiusCollision(radiusCollision_), orientation(orientation_) {
    pos[0] = pos_[0];
    pos[1] = pos_[1];
    pos[2] = pos_[2];
    Emap3(vel).setZero();
    Emap3(omega).setZero();
    sphNeighbor.reserve(10);
    return;
}

Sphere::~Sphere() noexcept {
    for (auto &l : sphLayer) {
        delete l.second;
        l.second = nullptr;
    }
    sphLayer.clear();
    sphNeighbor.clear();
}

Sphere::Sphere(const Sphere &other) noexcept {
    gid = other.gid;
    radius = other.radius;
    radiusCollision = other.radiusCollision;
    for (int i = 0; i < 3; i++) {
        pos[i] = other.pos[i];
        vel[i] = other.vel[i];
        omega[i] = other.omega[i];
    }
    orientation = other.orientation;

    // copy neighbors
    sphNeighbor = other.sphNeighbor;

    // copy layers
    for (auto &l : other.sphLayer) {
        sphLayer[l.first] = new SPHExp(*(l.second));
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
    sphNeighbor = std::move(other.sphNeighbor);
    sphLayer = std::move(other.sphLayer);
}

Sphere &Sphere::operator=(Sphere other) noexcept {
    using std::swap;
    swap(*this, other);
    return *this;
}

void Sphere::addLayer(const std::string &name, SPHExp::KIND kind, int order, const Equatn orientation) {
    // the orientation of the layer can be different from the orientation of the sphere
    sphLayer[name] = new SPHExp(kind, name, order, orientation);
}

void Sphere::dumpSphere() const {
    printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", gid, radius, radiusCollision, pos[0], pos[1], pos[2]);
    printf("vel %8f, %8f, %8f; omega %8f, %8f, %8f\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %8f, %8f, %8f, %8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
}

void Sphere::dumpNeighbor() const {
    for (auto &nb : sphNeighbor) {
        printf("gid %8d, r %8f, rCol %8f, pos %8f, %8f, %8f\n", nb.gid, nb.radius, nb.radiusCollision, nb.pos[0],
               nb.pos[1], nb.pos[2]);
        printf("Orientation: %8f, %8f, %8f, %8f; posRelative %8f, %8f, %8f\n", nb.orientation.w(), nb.orientation.x(),
               nb.orientation.y(), nb.orientation.z(), nb.posRelative[0], nb.posRelative[1], nb.posRelative[2]);
    }
}

void Sphere::dumpLayer(const std::string &name) const { sphLayer.find(name)->second->dumpSpectralValues(); }

void Sphere::Pack(std::vector<char> &buff) const {
    Buffer buffer(buff);
    // head
    buffer.pack(std::string("SPHERE"));
    // fixed size data
    buffer.pack(gid);                                                 // int gid = INVALID;
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
        const Equatn &orientation = layer.second->orientation;
        buffer.pack(std::array<double, 4>{orientation.w(), orientation.x(), orientation.y(),
                                          orientation.z()}); // Equatn orientation;
        buffer.pack(layer.second->spectralCoeff);
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
    sphNeighbor.clear();

    int nLayer = 0;
    buffer.unpack(nLayer, buff); // number of layers
    std::string strbuf2;
    sphLayer.clear();
    for (int i = 0; i < nLayer; i++) {
        buffer.unpack(strbuf, buff); // the name
        // data
        int kind, order;
        buffer.unpack(kind, buff);    //    const KIND kind;
        buffer.unpack(order, buff);   // const int order;        // the order of expansion p
        buffer.unpack(strbuf2, buff); //    const std::string name; // the name of this quantity
        buffer.unpack(array4, buff);  // Equatn orientation;
        Equatn orientTemp;
        orientTemp.w() = array4[0];
        orientTemp.x() = array4[1];
        orientTemp.y() = array4[2];
        orientTemp.z() = array4[3];
        // allocate
        sphLayer[strbuf] = new SPHExp(static_cast<SPHExp::KIND>(kind), strbuf2, order, orientTemp);
        // unpack the sph values
        buffer.unpack(sphLayer[strbuf]->spectralCoeff, buff);
    }
}
