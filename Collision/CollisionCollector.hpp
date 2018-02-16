#ifndef COLLISIONCOLLECTOR_HPP
#define COLLISIONCOLLECTOR_HPP

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

#include <omp.h>

#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"

struct CollisionBlock { // the information for each collision
  public:
    double phi0;  // constraint value
    double gamma; // force magnitude , could be an initial guess
    int gidI, gidJ;
    int globalIndexI, globalIndexJ;
    Evec3 normI, normJ; // norm vector for each particle. gvecJ = - gvecI
    Evec3 posI, posJ;   // the collision position on I and J. useless for spheres.

    CollisionBlock() : gidI(0), gidJ(0), globalIndexI(0), globalIndexJ(0), phi0(0), gamma(0) {
        // default constructor
        normI.setZero();
        normJ.setZero();
        posI.setZero();
        posJ.setZero();
    }

    CollisionBlock(double phi0_, double gamma_, int gidI_, int gidJ_, int globalIndexI_, int globalIndexJ_,
                   const Evec3 &normI_, const Evec3 &normJ_, const Evec3 &posI_, const Evec3 &posJ_)
        : phi0(phi0_), gamma(gamma_), gidI(gidI_), gidJ(gidJ_), globalIndexI(globalIndexI_),
          globalIndexJ(globalIndexJ_), normI(normI_), normJ(normJ_), posI(posI_), posJ(posJ_) {}
};

using CollisionBlockQue = std::vector<CollisionBlock>; // can be changed to other containers, e.g., deque
using CollisionBlockPool = std::vector<CollisionBlockQue>;

// process each collision i,j and record them to ColBlocks
// interface operator(obj i, obj j)
class CollisionCollector {
  public:
    std::shared_ptr<CollisionBlockPool> collisionPoolPtr;

    // constructor
    CollisionCollector() {
        const int totalThreads = omp_get_max_threads();
        collisionPoolPtr = std::make_shared<CollisionBlockPool>();
        collisionPoolPtr->resize(totalThreads);
        std::cout << "stress recoder size:" << collisionPoolPtr->size() << std::endl;
    }

    // copy constructor
    CollisionCollector(const CollisionCollector &obj) = default;
    CollisionCollector(CollisionCollector &&obj) = default;
    CollisionCollector &operator=(const CollisionCollector &obj) = default;
    CollisionCollector &operator=(CollisionCollector &&obj) = default;
    ~CollisionCollector() = default;

    void clear() {
        for (int i = 0; i < collisionPoolPtr->size(); i++) {
            (*collisionPoolPtr)[i].clear();
        }
    }

    template <class Trg, class Src>
    void operator(Trg &trg, Src &src) {
        const int threadId = omp_get_thread_num();
        auto &colque = (*collisionPoolPtr)[threadId];
        // construct a collision block to threadId
        CollisionBlock block;
        bool collide = trg.collide(src, block);
        if (collide) {
            colque.push_back(block);
        }
    }
};

#endif