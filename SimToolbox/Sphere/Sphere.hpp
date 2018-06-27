#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "Shexp.hpp"
#include "Util/Buffer.hpp"
#include "Util/GeoCommon.h"
#include "Util/IOHelper.hpp"

class Sphere {
  public:
    int gid = GEO_INVALID_INDEX;
    int globalIndex = GEO_INVALID_INDEX;
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

    void addLayer(const std::string &name_, const Shexp::KIND &kind_, const int &order_ = 4, const double &radius_ = -1,
                  const Equatn orientation_ = Equatn::Identity());

    void dumpSphere() const;
    // void dumpNeighbor() const;
    void dumpLayer(const std::string &name) const;

    Shexp &getLayer(const std::string &name);
    const Shexp &getLayer(const std::string &name) const;

    // motion
    void stepEuler(double dt);

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusCollision * 1.2; }

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
