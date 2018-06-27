#ifndef SYLINDER_HPP
#define SYLINDER_HPP

#include <unordered_map>
#include <vector>

#include "Util/Buffer.hpp"
#include "Util/EigenDef.hpp"
#include "Util/IOHelper.hpp"

#ifndef OBJ_INVALID_INDEX
#define OBJ_INVALID_INDEX (-1)
#endif

class Sylinder {
  public:
    int gid = OBJ_INVALID_INDEX;
    int globalIndex = OBJ_INVALID_INDEX;
    double radius;
    double radiusCollision;
    double length;
    double lengthCollision;
    double radiusSearch;
    Evec3 pos;
    Evec3 vel;
    Evec3 omega;
    Equatn orientation;

    // these are not packed and transferred
    Evec3 velCol;
    Evec3 omegaCol;
    Evec3 velBrown;
    Evec3 omegaBrown;
    Evec3 velHydro;
    Evec3 omegaHydro;

    Sylinder() = default;
    ~Sylinder() = default;

    Sylinder(const int &gid_, const double &radius_, const double &radiusCollision_, const double &length_,
             const double &lengthCollision_, const Evec3 &pos_ = Evec3::Zero(),
             const Equatn &orientation_ = Equatn::Identity());

    Sylinder(const Sylinder &) = default;
    Sylinder &operator=(const Sylinder &) = default;

    Sylinder(Sylinder &&) = default;
    Sylinder &operator=(Sylinder &&) = default;

    void dumpSylinder() const;

    void clear();

    // motion
    void stepEuler(double dt);

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusCollision * 4; }

    void Pack(std::vector<char> &buff) const;

    void Unpack(const std::vector<char> &buff);

    static void writeVTP(const std::vector<Sylinder> &sylinder, const std::string &prefix, const std::string &postfix,
                         int rank);
    static void writePVTP(const std::string &prefix, const std::string &postfix, const int nProcs);
};

#endif