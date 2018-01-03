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

void pack(Buffer & buffer) {
    // head
    buffer.pack(std::string("SPHERE"));
    // fixed size data
    buffer.pack(gid);    // int gid = INVALID;
    buffer.pack(radius); // double radius;
    buffer.pack(radiusCollision); // double radiusCollision;
    buffer.pack(std::array<double,3>(pos[0],pos[1],pos[2])); // Evec3 pos;
    buffer.pack(std::array<double,3>(vel[0],vel[1],vel[2])); // Evec3 vel;
    buffer.pack(std::array<double,3>(omega[0],omega[1],omega[2])); // Evec3 omega;
    buffer.pack(std::array<double,4>(orientation.w(),orientation.x(),orientation.y(),orientation.z())); // Equatn orientation;
    // variable size data
    const int nLayer=sphLayer.size();
    buffer.pack(nLayer); // number of layers
    for(auto &layer:sphLayer){
        // name
        buffer.pack(layer.first);
        // data
        buffer.pack(std::static_cast<int>(layer.second.kind)); //    const KIND kind;
        buffer.pack(layer.second.order); // const int order;        // the order of expansion p
        buffer.pack(layer.second.name); //    const std::string name; // the name of this quantity
        const Equatn & orientation=layer.second.orientation;
        buffer.pack(std::array<double,4>(orientation.w(),orientation.x(),orientation.y(),orientation.z())); // Equatn orientation;
        buffer.pack(layer.second.spectralCoeff);
    }
}

void unpack(const Buffer & buffer) {
    // head
    std::string strbuf;
    buffer.unpack(strbuf);    assert(head_==std::string("SPHERE"));
    // fixed size data
    buffer.unpack(gid);    // int gid = INVALID;
    buffer.unpack(radius); // double radius;
    buffer.unpack(radiusCollision); // double radiusCollision;
    std::array<double,3> array3;
    buffer.unpack(array3); // Evec3 pos;
    pos[0]=array3[0];
    pos[1]=array3[1];
    pos[2]=array3[2];
    buffer.unpack(array3); // Evec3 vel;
    vel[0]=array3[0];
    vel[1]=array3[1];
    vel[2]=array3[2];
    buffer.unpack(array3); // Evec3 omega;
    omega[0]=array3[0];
    omega[1]=array3[1];
    omega[2]=array3[2];
    std::array<double,4> array4;
    buffer.unpack(array4); // Equatn orientation;
    orientation[0]=array4[0];
    orientation[1]=array4[1];
    orientation[2]=array4[2];
    orientation[3]=array4[3];

    // variable size data
    sphNeighbor.clear();

    int nLayer=0;
    buffer.unpack(nLayer); // number of layers
    std::string strbuf2;
    sphLayer.clear();
    for(int i=0;i<nLayer;i++){
        buffer.unpack(strbuf);  // the name
        // data
        int kind,order;
        buffer.unpack(kind); //    const KIND kind;
        buffer.unpack(order); // const int order;        // the order of expansion p
        buffer.unpack(strbuf2); //    const std::string name; // the name of this quantity
        buffer.unpack(array4); // Equatn orientation;
        Equatn orientTemp;
        orientTemp[0]=array4[0];
        orientTemp[1]=array4[1];
        orientTemp[2]=array4[2];
        orientTemp[3]=array4[3];
        // allocate
        sphLayer[strbuf]= new SPHExp(std::static_cast<SPHExp::KIND>(kind), strbuf2, order, orientTemp);
        // unpack the sph values
        buffer.unpack(sphLayer[strbuf]->spectralCoeff);
    }

}
