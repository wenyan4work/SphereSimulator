#include <cstdio>
#include "Sphere.hpp"

void swap(Sphere &sphA, Sphere &sphB) {
    using std::swap;
    swap(sphA.gid, sphB.gid);
    swap(sphA.radius, sphB.radius);
    swap(sphA.radiusCollision, sphB.radiusCollision);

    sphA.pos.swap(sphB.pos);
    sphA.vel.swap(sphB.vel);
    sphA.omega.swap(sphB.omega);

    Equatn otemp = sphB.orientation;
    sphB.orientation = sphA.orientation;
    sphA.orientation = otemp;

    // swap the containers
    sphA.sphLayer.swap(sphB.sphLayer);
    sphA.sphNeighbor.swap(sphB.sphNeighbor);

    return;
}

Sphere::Sphere(int gid_, double radius_, double radiusCollision_) noexcept : gid(gid_),
                                                                             radius(radius_),
                                                                             radiusCollision(radiusCollision_) {
    pos.setZero();
    vel.setZero();
    omega.setZero();
    orientation.setIdentity();
    sphNeighbor.reserve(50);
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
    pos = other.pos;
    vel = other.vel;
    omega = other.omega;
    orientation = other.orientation;

    // copy neighbors
    sphNeighbor = other.sphNeighbor;

    // copy layers
    for (auto &l : other.sphLayer) {
        // TODO:copy the SPHobject
        sphLayer[l.first] = new SPHExp(*(l.second));
    }
}

Sphere::Sphere(Sphere &&other) noexcept {
    gid = other.gid;
    radius = other.radius;
    radiusCollision = other.radiusCollision;
    pos = other.pos;
    vel = other.vel;
    omega = other.omega;
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
    printf("orient %8f,%8f,%8f,%8f\n", orientation.w(), orientation.x(), orientation.y(), orientation.z());
}

void Sphere::dumpLayer(const std::string &name) const {
    // TODO: write this dump function for debugging
    sphLayer.find(name)->second->dumpSpectralValues();
}

void pack(char *const ptr) {}

void unpack(const char *const ptr) {}
