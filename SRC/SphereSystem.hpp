#ifndef SPHERESYSTEM_H_
#define SPHERESYSTEM_H_

#include "Config.h"
#include "MPI/InteractionManager.hpp"
#include "Sphere/CollisionSphere.hpp"
#include "Sphere/Sphere.hpp"

#include "STKFMM/STKFMM.h"

#include "SimToolbox/Collision/CollisionSolver.hpp"
#include "SimToolbox/Trilinos/TpetraUtil.hpp"
#include "SimToolbox/Util/IOHelper.hpp"
#include "SimToolbox/Util/TRngPool.hpp"

#include <algorithm>
#include <deque>
#include <memory>
#include <set>

#define COLBUF 0.3 // possible collision within (1+0.3)*(radiusI+radiusJ)

class SphereSystem {

  public:
    const Config runConfig;
    int snapID; // id of snap shot
    int stepCount;

    SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv);
    ~SphereSystem() = default;
    // forbid copy
    SphereSystem(const SphereSystem &) = delete;
    SphereSystem &operator=(const SphereSystem &) = delete;

    void moveEuler(); // Euler step update position and orientation, with given velocity
    void resolveCollision(bool manybody, double buffer = 0); // resolve collision
    void step();
    void output();    // output
    void partition(); // loadbalancing of spheres

  private:
    std::vector<Sphere> sphere; // spheres

    std::shared_ptr<TRngPool> rngPoolPtr; // thread safe rng

    std::shared_ptr<InteractionManager<double, 3, Sphere, Sphere>> interactManagerPtr; // near neighbor

    std::shared_ptr<CollisionSolver> collisionSolverPtr; // collision
    std::shared_ptr<CollisionCollector> collisionCollectorPtr;

    std::shared_ptr<stkfmm::STKFMM> fmmPtr;

    // MPI stuff
    Teuchos::RCP<const TCOMM> commRcp;       // mpi communicator
    Teuchos::RCP<TMAP> sphereMapRcp;         // 1 dof per sphere
    Teuchos::RCP<TMAP> sphereMobilityMapRcp; // 6 dof per sphere

    void setInitial(const std::string &initPosFile); // initial configuration of spheres
    void prepareTimestep();                          // prepare one timestep
    void statistics();                               // statistics

    bool readXYZ(const std::string &filename);

    // IO
    void writeVTK(const std::string &baseFolder);
    void writeXYZ(const std::string &baseFolder);
    std::vector<IOHelper::FieldVTU> dataFieldVTU;

    // TODO: implement these two for restart a simulation exactly
    void writeSerialized();
    void readSerialized();

    Teuchos::RCP<TOP> getMobOperator(bool manybody, std::string name);
    Teuchos::RCP<TV> getVelocityKnown(Teuchos::RCP<TOP> &mobilityOpRcp, Teuchos::RCP<TV> &forceRcp) const;
    Teuchos::RCP<TV> getVelocityBrown() const;
    Teuchos::RCP<TV> getForceKnown();

    void calcBoundaryCollision();
    void fitFMMBox();

    void writeBackVelocity(Teuchos::RCP<TV> &velocityRcp);

    void applyMonoLayer();
};

#endif
