/*
 * LayerSphereSyetem.h
 *
 *  Created on: Jan 2, 2018
 *      Author: wyan
 */

#ifndef HYDROSPHERE_H_
#define HYDROSPHERE_H_

#include <vector>
#include "mpi.h"

#include "EigenDef.hpp"
#include "FMMWrapper.h"
#include "TpetraUtil.hpp"
#include "ZDD.hpp"
#include "Sphere.hpp"

namespace LayerSphere {

// non dimensionalized
class LayerSphereSystem {
  public:
    enum SolverType {
        EXPLICIT,
        IMPLICIT
    };
    enum PrecType {
        RELAXATION,
        ILUT,
        KINV
    };

  private:
    // physics domain information
    Evec3 boxLow;
    Evec3 boxHigh;

    // solver information
    double iteres;
    int totalChebN;
    int precChoice; // choice of preconditioner

    // internal sphere data structure
    std::deque<Sphere> sphere;
    std::vector<SpherePacked> sendBuff;
    std::vector<SpherePacked> recvBuff;

    // FMM data structure
    std::vector<double> src_coord;
    std::vector<double> src_value;
    std::vector<double> trg_coord;
    std::vector<double> trg_value;
    FMM_Wrapper myFMM;

    // Trilinos data structure
    const Teuchos::RCP<const TCOMM> commRcp;
    ZDD<int> sphereGidFindDD;

    // MAP
    Teuchos::RCP<const TMAP> fdistMapRcp;     // contiguous fdist map across all nodes
    Teuchos::RCP<const TMAP> fdistFullMapRcp; // fully repeated fdist map across all nodes.
    Teuchos::RCP<const TMAP> mobilityMapRcp;  // a contiguous map for 6 entry per sphere

    Teuchos::RCP<TCMAT> approxImpSolverMatrixRcp; // the self+neighbor approximation to the full dense matrix
    Teuchos::RCP<TCMAT> flowRegMatrixRcp;         // regularize matrix with flow
    Teuchos::RCP<TV> fdistRcp;                    // force distribution
    Teuchos::RCP<TV> fdistTempSpaceRcp;           // force distribution, temp space durin calculation
    Teuchos::RCP<TV> bRcp;                        // the right side of implicit

    // MobilityMatrix
    // Teuchos::RCP<TCMAT> mobilityMatrixRcp;

    // known velocity and omega due to background flow, 6 unknowns per sphere
    Teuchos::RCP<TV> mobilityVelocityBackgroundFlowRcp;

    Teuchos::RCP<TOP> mobilityOperatorRcp;

    void exchangeSphere(); // move the old data according to the input sphere data.

    void locateSphere();

    void sortSphere(); // sort sphere according to sphereIO, must exchange before calling this.

    void packRigidSphere(const RigidSphere &, SpherePacked &);

    void unpackRigidSphere(const SpherePacked &, RigidSphere &);

    void syncSphereIn(); // copy sphereIO data to sphere, before solution

    void syncSphereInForceTorque(); // copy sphereIO force and torque to sphere

    void syncSphereOut(); // copy sphere data to sphereIO, after solution

    void updateNeighborBlocks(bool updateNbBlockMat);

    void updatefdistIndex();

    void updateFMMTree();

    void updateFdistWithBackgroundflow();

    void setFdistInitialGuess(Teuchos::RCP<TV> &); // setup initial guess with sphere.fdist, fdistlast

    void calculateSelfBlock();

    void calculateMobilityMatrix();

    // void calculateSelfBlockSingleSphere(RigidSphere &); // update self block of the solution matrix

    // void calculateMobilityMatrixSingleSphere(RigidSphere &); // calculate the mobility matrix of a single sphere

    // void calculateNeighborBlockSingleSphere(RigidSphere &, bool); // calculate the neighbor blocks of the iterative
    // matrix

    void constructDiscretization(); // refine mesh on each sphere with known meshDelta

    // void assembleMobilityDiagonalMatrix(const double, const double, const double);
    // assemble the block digonal matrix without neighbor blocks

    void assembleFlowRegMatrix(double droptol = 0);

    void assembleApproxImpSolverMatrix(double droptol = 0);

    void calcSphereMotionWithFdist(SolverType);

    void saveFdist(); // save fdist to fdistlast

    void calcBackgroundFlow(); // update background flow with fdist, use fdistRcp and fdistTempSpaceRcp as temp space

    void reSampleEvecXYZ(Evec &, const int, const ChebNodal *);

    ChebNodal *getChebIntegrator(const double &, const double &);

  public:
    const int myRank;
    const int nProcs;
    const double stkReg; // the regularization parameter for stokes kernel
    const double KReg;   // the regularization parameter for nonlocal K kernel

    // the input data to the system, specifying the current location, driving force, neighbor info, etc
    std::deque<RigidSphereIO> sphereIO;
    std::deque<pointFlow> pointFlowList;

    // max N, periodic type, fmm parameter, iterative solver parameter
    HydroSphereSystem(const int, const FMM_Wrapper::PAXIS, const int, const double, const Evec3 &, const Evec3 &,
                      const double, const double);

    ~HydroSphereSystem() = default;

    void addSphere(const RigidSphereIO &);

    void getReadyForSolve(SolverType); // must be called before calling the solve() functions

    void applyImpSolverMatrixFdist(const TMV &in, TMV &out, double alpha, double beta);

    void setPrecType(PrecType);

    void dumpSphere(int) const;

    void dumpSphereIO(int) const;

    void solveFullImplicit(); // implicitly solve background flow with density

    void solveExplicitWithBackgroundFlow(); // explicitly solve with known background flow

    void solveMobilityVelocityBackgroundFlow(); // for use with LCP collision code
    void updateFdistWithForceTorque();          // for use with LCP collision code

    void calcFlowOnGridWithFdist(double Lscale, double Tscale, double dx, std::string filename, double xShift = 0,
                                 double yShift = 0, double zShift = 0);

    Teuchos::RCP<TOP> &getImplicitMobilityOperator(const double, const double, const double); // get the mobility matrix

    Teuchos::RCP<TV> &getMobilityVelocityBackgroundFlow(const double, const double, const double);

    Teuchos::RCP<const TMAP> &getFdistMap();

    Teuchos::RCP<const TMAP> &getMobilityMap();
};
}

#endif /* HYDROFIBER_HYDRORIGID_H_ */
