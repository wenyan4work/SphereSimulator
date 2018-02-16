#ifndef COLLISIONSPHERE_HPP
#define COLLISIONSPHERE_HPP

#include "CollisionCollector.hpp"
#include "Sphere/Sphere.hpp"
#include "Util/Buffer.hpp"
#include "Util/EigenDef.hpp"

#define INVALID -1
#define COLBUF 0.3 // record collision block for sep < COLBUF*(rI+RJ)

class CollisionSphere {
    // this is a POD type, used to generate collision blocks
  public:
    int gid = INVALID;
    int globalIndex = INVALID;
    double radiusCollision;
    double pos[3];

    void CopyFromFull(const Sphere &s) {
        gid = s.gid;
        globalIndex = s.globalIndex;
        radiusCollision = s.radiusCollision;
        pos[0] = s.pos[0];
        pos[1] = s.pos[1];
        pos[2] = s.pos[2];
    }

    // necessary interface for Near Interaction
    const double *Coord() const { return pos; }
    double Rad() const { return radiusCollision * 2; }

    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        mybuff.pack(gid);
        mybuff.pack(globalIndex);
        mybuff.pack(radiusCollision);
        mybuff.pack(pos[0]);
        mybuff.pack(pos[1]);
        mybuff.pack(pos[2]);
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        mybuff.unpack(gid, buff);
        mybuff.unpack(globalIndex, buff);
        mybuff.unpack(radiusCollision, buff);
        mybuff.unpack(pos[0], buff);
        mybuff.unpack(pos[1], buff);
        mybuff.unpack(pos[2], buff);
    }

    inline bool collide(const CollisionSphere &sphereJ, CollisionBlock &block) {
        if (gid == sphereJ.gid) {
            // no self collision
            return false;
        }
        EAvec3 posI(pos[0], pos[1], pos[2]);
        EAvec3 posJ(sphereJ.pos[0], sphereJ.pos[1], sphereJ.pos[2]);
        EAvec3 rIJ = posJ - posI;
        const double rIJNorm = rIJ.norm();
        const double sep = rIJNorm - radiusCollision - sphereJ.radiusCollision;

        // save collision block
        // save only block gidI < gidJ
        if (sep < COLBUF * (radiusCollision + sphereJ.radiusCollision) && gid < sphereJ.gid) {
            // collision
            block.normI = rIJ * (-1.0 / rIJNorm);
            block.normJ = -block.normI;
            block.phi0 = sep;
            block.gidI = gid;
            block.gidJ = sphereJ.gid;
            block.globalIndexI = globalIndex;
            block.globalIndexJ = sphereJ.globalIndex;
            block.posI.setZero();
            block.posJ.setZero();
            block.gamma = sep < 0 ? -sep : 0; // a crude initial guess
            return true;
        } else {
            // no collision
            return false;
        }
    }

    friend void swap(CollisionSphere &, CollisionSphere &);
};

void swap(CollisionSphere &A, CollisionSphere &B) {
    using std::swap;
    swap(A.gid, B.gid);
    swap(A.globalIndex, B.globalIndex);
    swap(A.radiusCollision, B.radiusCollision);
    swap(A.pos[0], B.pos[0]);
    swap(A.pos[1], B.pos[1]);
    swap(A.pos[2], B.pos[2]);
}

#endif