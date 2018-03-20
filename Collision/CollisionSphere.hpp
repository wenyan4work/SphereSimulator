#ifndef COLLISIONSPHERE_HPP
#define COLLISIONSPHERE_HPP

#include "CollisionCollector.hpp"
#include "Sphere/Sphere.hpp"
#include "Util/Buffer.hpp"
#include "Util/EigenDef.hpp"

#define INVALID -1
#define COLBUF 0.3 // record collision block for sep < COLBUF*(rI+RJ)

class CollisionSphere {
  public:
    int gid = INVALID;
    int globalIndex = INVALID;
    double radiusCollision;
    Evec3 pos;

    void CopyFromFull(const Sphere &s) {
        gid = s.gid;
        globalIndex = s.globalIndex;
        radiusCollision = s.radiusCollision;
        pos = s.pos;
    }

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusCollision * (1 + COLBUF * 2); }

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
        if (gid >= sphereJ.gid) {
            // no self collision
            // do not record gid > J.gid
            return false;
        }
        EAvec3 posI(pos[0], pos[1], pos[2]);
        EAvec3 posJ(sphereJ.pos[0], sphereJ.pos[1], sphereJ.pos[2]);
        EAvec3 rIJ = posJ - posI;
        const double rIJNorm = rIJ.norm();
        const double sep = rIJNorm - radiusCollision - sphereJ.radiusCollision;

        // save collision block
        // save only block gidI < gidJ
        if (sep < COLBUF * radiusCollision) {
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
};

#endif