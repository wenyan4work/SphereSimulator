#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "Shexp.hpp"
#include "Util/Buffer.hpp"
#include "Util/IOHelper.hpp"

#define INVALID -1

// class NeighborSphere {
//   public:
//     int gid = INVALID;
//     double radius;
//     double radiusCollision;
//     Evec3 pos;
//     Evec3 posRelative;
//     Equatn orientation;
//     NeighborSphere() = default;
//     NeighborSphere(const int &gid_, const Evec3 &posRelative_) noexcept : gid(gid_), posRelative(posRelative_) {}
// };

class Sphere {
  public:
    int gid = INVALID;
    int globalIndex = INVALID;
    double radius;
    double radiusCollision;
    Evec3 pos;
    Evec3 vel;
    Evec3 omega;
    Equatn orientation;

    // std::vector<NeighborSphere> sphNeighbor;
    std::unordered_map<std::string, Shexp *> sphLayer;

    Sphere() = default;
    Sphere(const int &gid_, const double &radius_, const double &radiusCollision_, const Evec3 &pos_ = Evec3::Zero(),
           const Equatn &orientation_ = Equatn::Identity()) noexcept;
    ~Sphere() noexcept;

    Sphere(const Sphere &) noexcept;
    Sphere(Sphere &&) noexcept;

    Sphere &operator=(Sphere) noexcept;

    void addLayer(const std::string &, Shexp::KIND, int order = 4, const Equatn orientation = Equatn::Identity());

    void dumpSphere() const;
    // void dumpNeighbor() const;
    void dumpLayer(const std::string &name) const;

    Shexp &getLayer(const std::string &name);
    const Shexp &getLayer(const std::string &name) const;

    // motion
    void stepEuler(double dt);

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusCollision * 2; }

    void Pack(std::vector<char> &buff) const;

    void Unpack(const std::vector<char> &buff);

    friend void swap(Sphere &, Sphere &);

    static void writeVTP(const std::vector<Sphere> &sphere, const std::string &prefix, const std::string &postfix,
                         int rank);
    static void writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs);

    static void writeVTU(const std::vector<Sphere> &sphere, const std::vector<IOHelper::FieldVTU> &dataFields,
                         const std::string &prefix, const std::string &postfix, int rank);
    static void writePVTU(const std::vector<IOHelper::FieldVTU> &dataFields, const std::string &prefix,
                          const std::string &postfix, const int nProcs);
};

#endif // SPHERE_HPP
