#ifndef EQUATNHELPER_HPP
#define EQUATNHELPER_HPP

// add a few functions to eigen quaternion

#include <cmath>

#include "EigenDef.hpp"

// WARNING: use w,x,y,z of a quaternion, do not use [0,1,2,3], since the mapping order is different in eigen and in most
// papers. eigen: w->[3], for scalar. x,y,z->[0,1,2] for vector

// https://eigen.tuxfamily.org/dox/group__TutorialGeometry.html

// Equatn: This class represents a quaternion $ w+xi+yj+zk $

// eigen builtin routine : FromTwoVectors(a,b)
// the built rotation represent a rotation sending the line of direction a to the line of direction b, both lines
// passing through the origin. Note that the two input vectors do not have to be normalized, and do not need to have the
// same norm.

// Concatenation of two transformations
// gen1 * gen2; Apply the transformation to a vector
// vec2 = gen1 * vec1; Get the inverse of the transformation
// gen2 = gen1.inverse(); Spherical interpolation  (Rotation2D and Quaternion only)

class EquatnHelper {
  public:
    static void setUnitRandomEquatn(Equatn &q, const double &u1, const double &u2, const double &u3) {
        // a random unit quaternion following a uniform distribution law on SO(3)
        // from three U[0,1] random numbers
        constexpr double pi = 3.14159265358979323846;
        const double a = sqrt(1 - u1);
        const double b = sqrt(u1);
        const double su2 = sin(2 * pi * u2);
        const double cu2 = cos(2 * pi * u2);
        const double su3 = sin(2 * pi * u3);
        const double cu3 = cos(2 * pi * u3);
        q.w() = a * su2;
        q.x() = a * cu2;
        q.y() = b * su3;
        q.z() = b * cu3;
    }

    static void rotateEquatn(Equatn &q, const Evec3 &omega, const double &dt) {
        // Delong, JCP, 2015, Appendix A eq1, not linearized
        const double w = omega.norm();
        if (w < std::numeric_limits<float>::epsilon()) {
            return;
        }
        const double winv = 1 / w;
        const double sw = sin(w * dt / 2);
        const double cw = cos(w * dt / 2);
        const double s = q.w();
        const EAvec3 p(q.x(), q.y(), q.z());
        const EAvec3 xyz = s * sw * omega * winv + cw * p + sw * winv * (omega.cross(p));
        q.w() = s * cw - (p.dot(omega)) * sw * winv;
        q.x() = xyz.x();
        q.y() = xyz.y();
        q.z() = xyz.z();
        q.normalize();
    }

    static void getCrossProductMatrix(const Evec3 &p, Emat3 &P) {
        P(0, 0) = 0;
        P(0, 1) = -p[2];
        P(0, 2) = p[1];
        P(1, 0) = p[2];
        P(1, 1) = 0;
        P(1, 2) = -p[0];
        P(2, 0) = -p[1];
        P(2, 1) = p[0];
        P(2, 2) = 0;
        return;
    }

    static void getPsiMatFromEquatn(const Equatn &q, EmatPsi &psi) {
        const double s = q.w();
        const EAvec3 p(q.x(), q.y(), q.z());
        psi.block<1, 3>(0, 0) = 0.5 * (-p.transpose());
        psi(1, 0) = 0.5 * s;
        psi(1, 1) = 0.5 * p[2];
        psi(1, 2) = -0.5 * p[1];
        psi(2, 0) = -0.5 * p[2];
        psi(2, 1) = 0.5 * s;
        psi(2, 2) = 0.5 * p[0];
        psi(3, 0) = 0.5 * p[1];
        psi(3, 1) = -0.5 * p[0];
        psi(3, 2) = 0.5 * s;
    }
};

#endif
