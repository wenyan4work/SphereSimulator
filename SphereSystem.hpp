#ifndef SPHERESYSTEM_H_
#define SPHERESYSTEM_H_

#include <algorithm>
#include <deque>
#include <memory>
#include <set>

#include "Collision/CollisionSphere.hpp"
#include "MPI/InteractionManager.hpp"
#include "Sphere/Sphere.hpp"
#include "Util/TRngPool.hpp"
#include "STKFMM/STKFMM.h"

#include "config.h"

#define COLBUF 0.3 // possible collision within (1+0.3)*(radiusI+radiusJ)

class SphereSystem {

  public:
    double sysTime;
    double snapTime;
    int id_snap;

    void stepEuler(); // Euler step forward

    void output(); // output

    SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv);
    ~SphereSystem();
    // forbid copy
    SphereSystem(const SphereSystem &) = delete;
    SphereSystem &operator=(const SphereSystem &) = delete;

  private:
    const Config runConfig;
    // spheres
    std::vector<Sphere> sphere;
    // thread safe rng
    std::shared_ptr<TRngPool> rngPoolPtr;
    // MPI stuff
    Teuchos::RCP<TCOMM> commRcp;
    Teuchos::RCP<TMAP> sphereMapRcp; // 1 dof per sphere

    void setInitial(const std::string &initPos); // initial configuration of spheres
    void prepareTimestep();                      // prepare one timestep
    void statistics();                           // statistics
    void partition();                            // loadbalancing of objects

    // IO
    void writeVTK();
    void readVTK();

    Teuchos::RCP<TOP> getMobOperator(bool manybody);
    Teuchos::RCP<TOP> getVelocityKnown();
};

//#pragma GCC pop_options

#endif
