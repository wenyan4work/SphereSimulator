#ifndef COLLISIONSPHERE_HPP
#define COLLISIONSPHERE_HPP

#include "CollisionCollector.hpp"
#include "DCPQuery.hpp"
#include "Sylinder/Sylinder.hpp"
#include "Util/Buffer.hpp"
#include "Util/EigenDef.hpp"

constexpr double COLBUF = 0.5; // record collision block for sep < COLBUF*(rI+rJ)

class CollisionSylinder {
    // this is a POD type, used to generate collision blocks
  public:
    int gid = INVALID;
    int globalIndex = INVALID;
    double radiusCollision;
    double lengthCollision;
    double radiusSearch;
    Evec3 pos;
    Evec3 direction;

    DCPQuery<3, double, EAvec3> DistSegSeg3; // not transferred over mpi

    void CopyFromFull(const Sylinder &s) {
        gid = s.gid;
        globalIndex = s.globalIndex;
        radiusCollision = s.radiusCollision;
        lengthCollision = s.lengthCollision;
        radiusSearch = s.radiusSearch;
        pos = s.pos;
        direction = s.orientation * Evec3(0, 0, 1);
    }

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusSearch; }

    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        mybuff.pack(gid);
        mybuff.pack(globalIndex);
        mybuff.pack(radiusCollision);
        mybuff.pack(lengthCollision);
        mybuff.pack(radiusSearch);
        mybuff.pack(pos[0]);
        mybuff.pack(pos[1]);
        mybuff.pack(pos[2]);
        mybuff.pack(direction[0]);
        mybuff.pack(direction[1]);
        mybuff.pack(direction[2]);
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        mybuff.unpack(gid, buff);
        mybuff.unpack(globalIndex, buff);
        mybuff.unpack(radiusCollision, buff);
        mybuff.unpack(lengthCollision, buff);
        mybuff.unpack(radiusSearch, buff);
        mybuff.unpack(pos[0], buff);
        mybuff.unpack(pos[1], buff);
        mybuff.unpack(pos[2], buff);
        mybuff.unpack(direction[0], buff);
        mybuff.unpack(direction[1], buff);
        mybuff.unpack(direction[2], buff);
    }

    inline bool collide(const CollisionSylinder &sJ, CollisionBlock &block) {
        if (gid >= sJ.gid) {
            // no self collisio
            // do not record gid > J.gidn
            return false;
        }
        const auto &sI = *this;
        EAvec3 posI = sI.pos;
        EAvec3 posJ = sJ.pos;

        const EAvec3 Pm = posI - (0.5 * sI.lengthCollision) * sI.direction;
        const EAvec3 Pp = posI + (0.5 * sI.lengthCollision) * sI.direction;

        const EAvec3 Qm = posJ - (0.5 * sJ.lengthCollision) * sJ.direction;
        const EAvec3 Qp = posJ + (0.5 * sJ.lengthCollision) * sJ.direction;

        // location
        EAvec3 Ploc = Evec3::Zero(), Qloc = Evec3::Zero();
        double s = 0, t = 0;

        // check collision for line segment I and J
        const double distMin = DistSegSeg3(Pm, Pp, Qm, Qp, Ploc, Qloc, s, t);
        const double sep = distMin - (sI.radiusCollision + sJ.radiusCollision);

        // save collision block
        // save only block gidI < gidJ
        if (sep < COLBUF * (radiusCollision + sJ.radiusCollision)) {
            // collision
            block.normI = (Ploc - Qloc).normalized();
            block.normJ = -block.normI;
            block.phi0 = sep;
            block.gidI = gid;
            block.gidJ = sJ.gid;
            block.globalIndexI = globalIndex;
            block.globalIndexJ = sJ.globalIndex;
            block.posI = Ploc - posI; // collision location relative to geometric center
            block.posJ = Qloc - posJ;
            block.gamma = sep < 0 ? -sep : 0; // a crude initial guess
            return true;
        } else {
            // no collision
            return false;
        }
    }
};

#endif