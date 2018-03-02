#ifndef SPHERESYSTEM_H_
#define SPHERESYSTEM_H_

#include <algorithm>
#include <deque>
#include <memory>
#include <set>

#include "Collision/CollisionSolver.hpp"
#include "Collision/CollisionSphere.hpp"
#include "MPI/InteractionManager.hpp"
// #include "STKFMM/STKFMM.h"
#include "Sphere/Sphere.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/TRngPool.hpp"

#include "Config.h"

#define COLBUF 0.3 // possible collision within (1+0.3)*(radiusI+radiusJ)

class SphereSystem {

  public:
    const Config runConfig;
    double sysTime;
    double snapTime;
    int snapID; // id of snap shot

    SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv);
    ~SphereSystem() = default;
    // forbid copy
    SphereSystem(const SphereSystem &) = delete;
    SphereSystem &operator=(const SphereSystem &) = delete;

    void moveEuler(Teuchos::RCP<TV> &velocityRcp); // Euler step update position and orientation, with given velocity
    void resolveCollision(bool manybody, double buffer = 0); // resolve collision
    void step();
    void output();    // output
    void partition(); // loadbalancing of spheres

  private:
    std::vector<Sphere> sphere;           // spheres
    std::shared_ptr<TRngPool> rngPoolPtr; // thread safe rng

    std::shared_ptr<InteractionManager<double, 3, Sphere, Sphere>> interactManagerPtr;
    std::shared_ptr<CollisionSolver> collisionSolverPtr;
    std::shared_ptr<CollisionCollector> collisionCollectorPtr;

    // MPI stuff
    Teuchos::RCP<const TCOMM> commRcp;             // mpi communicator
    Teuchos::RCP<TMAP> sphereMapRcp;         // 1 dof per sphere
    Teuchos::RCP<TMAP> sphereMobilityMapRcp; // 6 dof per sphere

    void setInitial(const std::string &initPosFile); // initial configuration of spheres
    void prepareTimestep();                          // prepare one timestep
    void statistics();                               // statistics

    bool readXYZ(const std::string &filename);

    // IO
    void writeVTK(const std::string &baseFolder);
    void writeXYZ(const std::string &baseFolder);

    // TODO: implement these two for restart a simulation exactly
    void writeSerialized();
    void readSerialized();

    Teuchos::RCP<TOP> getMobOperator(bool manybody) const;
    Teuchos::RCP<TV> getVelocityKnown(Teuchos::RCP<TOP> &mobilityOpRcp, Teuchos::RCP<TV> &forceRcp) const;
    Teuchos::RCP<TV> getVelocityBrown() const;
    Teuchos::RCP<TV> getForceKnown() const;
};

#endif
