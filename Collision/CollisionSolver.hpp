#ifndef COLLISIONSOLVER_HPP
#define COLLISIONSOLVER_HPP

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"

#include "CPSolver.hpp"
#include "CollisionCollector.hpp"

class CollisionSolver {

  public:
    // constructor
    CollisionSolver() = default;
    ~CollisionSolver() = default;

    // forbid copy
    CollisionSolver(const CollisionSolver &) = delete;
    CollisionSolver(CollisionSolver &&) = delete;
    const CollisionSolver &operator=(const CollisionSolver &) = delete;
    const CollisionSolver &operator=(CollisionSolver &&) = delete;

    // working functions
    void reset() {
        res = 1e-5;
        maxIte = 2000;
        newton = false;

        objMobMapRcp.reset(); // distributed map for obj mobility. 6 dof per obj
        forceColRcp.reset();    // force vec, 6 dof per obj
        velocityColRcp.reset(); // velocity vec, 6 dof per obj. velocity = mobity * forceCol

        gammaMapRcp.reset(); // distributed map for collision magnitude gamma
        gammaRcp.reset();    // the unknown

        phi0Rcp.reset();
        vnRcp.reset();
        bRcp.reset(); // the constant piece of LCP problem

        matMobilityRcp.reset(); // mobility operator, 6 dof per obj to 6 dof per obj
        matFcTransRcp.reset();  // FcTrans matrix, 6 dof per obj to gamma dof

        queueThreadIndex.clear();
    }

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
    Teuchos::RCP<TV> getVelocityCol() const { return velocityColRcp; }
    Teuchos::RCP<TV> getForceColMagnitude() const { return gammaRcp; }
    Teuchos::RCP<TV> getVelocityKnown() const { return vnRcp; }
    Teuchos::RCP<TV> getPhi0() const { return phi0Rcp; }

  private:
    double res;
    int maxIte;
    bool newton;

    // mobility
    Teuchos::RCP<TMAP> objMobMapRcp; // distributed map for obj mobility. 6 dof per obj
    Teuchos::RCP<TV> forceColRcp;    // force vec, 6 dof per obj
    Teuchos::RCP<TV> velocityColRcp; // velocity vec, 6 dof per obj. velocity = mobity * forceCol

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