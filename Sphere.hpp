#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "SPHExp.hpp"
#include <unordered_map>

class Sphere {
  public:
    int gid;
    double radius;
    double radiusCollision;
    Evec3 pos;
    Evec3 vel;
    Evec3 omega;
    Equatn orientation;

    std::unordered_map<std::string, SPHExp *> sphPotential;

    Sphere() noexcept;
    ~Sphere() noexcept;

    Sphere(const Sphere &) noexcept;
    Sphere(const Sphere &&) noexcept;

    Sphere &operator=(Sphere) noexcept;

    friend void swap(Sphere &, Sphere &);

  private:
};

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

    // swap the unordered map
    sphA.sphPotential.swap(sphB.sphPotential);

    return;
}

#endif // SPHERE_HPP
