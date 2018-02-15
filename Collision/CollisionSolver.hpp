#ifndef COLLISIONSOLVER_HPP
#define COLLISIONSOLVER_HPP

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

#include "Collision/CPSolver.hpp"
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

class CollisionSolver {

  public:
    // constructor
    CollisionSolver();
    ~CollisionSolver();

    // forbid copy
    CollisionSolver(const CollisionSolver &) = delete;
    CollisionSolver(CollisionSolver &&) = delete;
    const CollisionSolver &operator=(const CollisionSolver &) = delete;
    const CollisionSolver &operator=(CollisionSolver &&) = delete;

    // working functions
    void setControlLCP(double res_, int maxIte_, bool newton_) {
        res = res_;
        maxIte = maxIte_;
        newton = newton_; // run mmNewton following the GD to refine the solution to 100x smaller residue
    }

    void setup(CollisionBlockPool &collision_, Teuchos::RCP<TMAP> &objMobMapRcp_, double dt_, double bufferGap_ = 0);

    // for implicit solver:
    // for explicit solver: pass the background velocity
    void solveCollision(Teuchos::RCP<TOP> &matMobilityRcp_, Teuchos::RCP<TV> &velocityKnownRcp_);

    // TODO: implement this in the future
    void solveCollisionSplit(Teuchos::RCP<TOP> &matMobilityMajorRcp, Teuchos::RCP<TOP> &matMobilityOtherRcp);

    // return results
    Teuchos::RCP<TV> getForceCol() const { return forceColRcp; }
    Teuchos::RCP<TV> getVelocityCol() const { return velocityRcp; }
    Teuchos::RCP<TV> getGammaCol() const { return gammaRcp; }

  private:
    double res = 1e-5;
    int maxIte = 2000;
    bool newton = false;

    Teuchos::RCP<TCOMM> commRcp;

    // mobility
    Teuchos::RCP<TMAP> objMobMapRcp; // distributed map for obj mobility. 6 dof per obj
    Teuchos::RCP<TV> forceColRcp;    // force vec, 6 dof per obj
    Teuchos::RCP<TV> velocityRcp;    // velocity vec, 6 dof per obj. velocity = mobity * forceCol

    // unknown collision force magnitude
    Teuchos::RCP<TMAP> gammaMapRcp; // distributed map for collision magnitude gamma
    Teuchos::RCP<TV> gammaRcp;      // the unknown

    // known initial value of constraints
    Teuchos::RCP<TV> phi0Rcp;
    Teuchos::RCP<TV> vnRcp;
    Teuchos::RCP<TV> bRcp; // the constant piece of LCP problem

    // Mobility operator and FcTrans matrices
    Teuchos::RCP<TOP> matMobilityRcp;  // mobility operator, 6 dof per obj to 6 dof per obj
    Teuchos::RCP<TCMAT> matFcTransRcp; // FcTrans matrix, 6 dof per obj to gamma dof

    std::vector<int> queueThreadIndex;

    void dumpCollision(CollisionBlockPool &collision_) const;
    void setupCollisionBlockQueThreadIndex(CollisionBlockPool &collision_);

    // utilities for FcTrans
    void setupFcTrans(CollisionBlockPool &collision_);
    // current known value of constraints
    void setupPhi0Vec(CollisionBlockPool &collision_, double dt_, double bufferGap_);
    // initial guess of unknown gamma
    void setupGammaVec(CollisionBlockPool &collision_);
    // known velocity
    void setupVnVec(CollisionBlockPool &collision_, std::vector<double> &velocity_);

    void setupBVec();

    int getNumNegElements(Teuchos::RCP<TV> &vecRcp_) const;
};

#endif