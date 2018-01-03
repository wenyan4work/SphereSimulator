#ifndef SPHSPHERESYSTEM_H_
#define SPHSPHERESYETEM_H_

#include "mpi.h"
#include <vector>

#include "Buffer.hpp"
#include "EigenDef.hpp"
#include "FMMWrapper.h"
#include "Sphere.hpp"
#include "TpetraUtil.hpp"
#include "ZDD.hpp"

#define SPHEREBYTESIZEMAX 48000
// default to 48KB per sphere.
// affect the buffer size for mpi send/recv 

class SPHSphereSystem {
    // system class provide common operations to all operators

  private:
    // physics domain information
    Evec3 boxLow;
    Evec3 boxHigh;

    // internal sphere data structure
    std::deque<Sphere> sphere;

    // FMM data structure
    std::vector<double> srcCoord;
    std::vector<double> srcValue;
    std::vector<double> trgCoord;
    std::vector<double> trgValue;
    FMM_Wrapper myFMM;

    // Trilinos data structure
    const Teuchos::RCP<const TCOMM> commRcp;
    ZDD<int> sphereGidFindDD;

    // MAP
    Teuchos::RCP<const TMAP> mobilityMapRcp; // a contiguous map for 6 entry per sphere
    Teuchos::RCP<const TMAP> sphereMapRcp;   // a contiguous map for 6 entry per sphere

    // MobilityMatrix
    // Teuchos::RCP<TCMAT> mobilityMatrixRcp;


    void locateSphere();   // find the rank id for each sphere
    void exchangeSphere(); // move the old data according to the input sphere data.
    void updateNeighborSphere(); // update the neighbor list of each sphere
    void updateMap();  // update the TMAP
    void sortSphere(); // sort sphere according to sphereIO, must exchange before calling this.

    // void syncSphereIn();  // copy sphereIO data to sphere, before solution
    // void syncSphereOut(); // copy sphere data to sphereIO,  after solution

  public:
    const int myRank;
    const int nProcs;

    // the input data to the system, specifying the current location, driving force, neighbor info, etc
    std::deque<SphereIO> sphereIO;

    // SPH layer order, periodic type, fmm parameter, iterative solver parameter
    SPHSphereSystem(const Evec3 &boxLow_, const Evec3 &boxHigh_, const int &approxNumberSphere,
                    const int &multOrder = 8, const FMM_Wrapper::PAXIS &pbc = FMM_Wrapper::PAXIS::NONE);

    // forbid copy
    SPHSphereSystem(const SPHSphereSystem &) = delete;
    SPHSphereSystem(SPHSphereSystem &&) = delete;
    SPHSphereSystem &operator=(const SPHSphereSystem &) = delete;
    SPHSphereSystem &operator=(SPHSphereSystem &&) = delete;

    ~SPHSphereSystem() = default;

    void addSphere(const SphereIO &);

    void getReadyForOperator(); // must be called before calling the solve() functions

    void dumpSphere() const;
    void dumpSphereIO() const;
};

#endif /* HYDROFIBER_HYDRORIGID_H_ */
