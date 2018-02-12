#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <unordered_map>
#include <vector>

#include "Buffer.hpp"
#include "SPHExp.hpp"

#define INVALID -1

class NeighborSphereIO {
  public:
    int gid = INVALID;
    Evec3 posRelative;
    NeighborSphereIO(const int &gid_, const Evec3 &posRelative_) noexcept : gid(gid_), posRelative(posRelative_) {}
};

class SphereIO {
  public:
    int gid = INVALID;
    double radius;
    double radiusCollision;
    Evec3 pos;
    Evec3 vel;
    Evec3 omega;
    Equatn orientation;
    std::vector<NeighborSphereIO> sphNeighborIO;

    SphereIO(const int &gid_, const double &radius_, const double &radiusCollision_, const Evec3 &pos_ = Evec3::Zero(),
             const Equatn &orientation_ = Equatn::Identity()) noexcept
        : gid(gid_), radius(radius_), radiusCollision(radiusCollision_), pos(pos_), orientation(orientation_) {
        vel.setZero();
        omega.setZero();
        sphNeighborIO.clear();
        sphNeighborIO.reserve(10);
    }

    void dumpSphere() const;
};

class NeighborSphere {
  public:
    int gid = INVALID;
    double radius;
    double radiusCollision;
    Evec3 pos;
    Evec3 posRelative;
    Equatn orientation;
    NeighborSphere() noexcept = default;
    NeighborSphere(const int &gid_, const Evec3 &posRelative_) noexcept : gid(gid_), posRelative(posRelative_) {}
};

class Sphere {
  public:
    int gid = INVALID;
    double radius;
    double radiusCollision;
    Evec3 pos;
    Evec3 vel;
    Evec3 omega;
    Equatn orientation;

    std::vector<NeighborSphere> sphNeighbor;
    std::unordered_map<std::string, SPHExp *> sphLayer;

    Sphere() = default;
    Sphere(const SphereIO& sphere) ;
    Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const Evec3 &pos_ = Evec3::Zero(),
           const Equatn &orientation_ = Equatn::Identity()) noexcept;
    ~Sphere() noexcept;

    Sphere(const Sphere &) noexcept;
    Sphere(Sphere &&) noexcept;

    Sphere &operator=(Sphere) noexcept;

    void addLayer(const std::string &, SPHExp::KIND, int order = 4, const Equatn orientation = Equatn::Identity());

    void dumpSphere() const;
    void dumpNeighbor() const;
    void dumpLayer(const std::string &) const;

    void pack(Buffer &buffer);
    void unpack(const Buffer &buffer);

    friend void swap(Sphere &, Sphere &);
};

#endif // SPHERE_HPP
