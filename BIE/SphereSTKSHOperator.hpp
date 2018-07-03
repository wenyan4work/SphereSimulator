#ifndef SPHERESTKSHOPERATOR_HPP_
#define SPHERESTKSHOPERATOR_HPP_

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "MPI/InteractionManager.hpp"
#include "STKFMM/STKFMM.h"
#include "Sphere/Sphere.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"
#include "Util/GeoCommon.h"

using namespace stkfmm;

class NearEvalSH {
  public:
    // int gid = GEO_INVALID_INDEX;
    double radiusNear;
    Evec3 pos;
    // grid value and sh data
    Shexp sh;
    // spectral value
    std::vector<double> spectralCoeff;
    // grid
    // grid center offset not applied to this grid Coord,
    // to be compatible with PBC images
    std::vector<double> gridWeight;
    std::vector<double> gridCoord;
    std::vector<double> gridNorm;

    // necessary interface for Near Interaction
    const double *Coord() const { return pos.data(); }

    double Rad() const { return radiusNear; }

    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        // mybuff.pack(gid);
        mybuff.pack(radiusNear);
        mybuff.pack(pos[0]);
        mybuff.pack(pos[1]);
        mybuff.pack(pos[2]);
        // pack sh
        mybuff.pack(static_cast<int>(sh.kind)); //    const KIND kind;
        mybuff.pack(sh.order);                  // const int order;        // the order of expansion p
        mybuff.pack(sh.name);
        mybuff.pack(sh.radius);
        const Equatn &orientation = sh.orientation;
        mybuff.pack(std::array<double, 4>{orientation.w(), orientation.x(), orientation.y(),
                                          orientation.z()}); // Equatn orientation;
        mybuff.pack(sh.gridValue);
        // pack other data
        mybuff.pack(spectralCoeff);
        mybuff.pack(gridWeight);
        mybuff.pack(gridCoord);
        mybuff.pack(gridNorm);
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        // mybuff.unpack(gid, buff);
        mybuff.unpack(radiusNear, buff);
        mybuff.unpack(pos[0], buff);
        mybuff.unpack(pos[1], buff);
        mybuff.unpack(pos[2], buff);
        int temp;
        mybuff.unpack(temp, buff);
        sh.kind = static_cast<Shexp::KIND>(temp);
        mybuff.unpack(sh.order, buff);
        mybuff.unpack(sh.name, buff);
        mybuff.unpack(sh.radius, buff);
        std::array<double, 4> orientation;
        mybuff.unpack(orientation, buff);
        sh.orientation.w() = orientation[0];
        sh.orientation.x() = orientation[1];
        sh.orientation.y() = orientation[2];
        sh.orientation.z() = orientation[3];
        mybuff.unpack(sh.gridValue, buff);
        mybuff.unpack(spectralCoeff, buff);
        mybuff.unpack(gridWeight, buff);
        mybuff.unpack(gridCoord, buff);
        mybuff.unpack(gridNorm, buff);
    }
};

class NearEvaluator {
  private:
    static std::shared_ptr<STKFMM> fmmPtr;

  public:
    // TODO: Stokes DL
    double cSL, cTrac;

    NearEvaluator(const double sl, const double trac) {
        cSL = sl;
        cTrac = trac;
    }

    void operator()(NearEvalSH &trg, NearEvalSH &src);
};

class SphereSTKSHOperator : public TOP {
  private:
    const double cId;
    const double cSL;
    const double cTrac;
    const double cLOP;

    const std::vector<Sphere> *const spherePtr; // read only
    const std::string name;                     // which sphere harmonics to work on

    // temporary data
    mutable std::vector<Shexp> sh;
    mutable std::vector<double> pointValues;      // a temporary space to store input vector
    mutable std::vector<double> pointValuesApply; // a temporary space to accumulate operator results

    // temporary data for sh coeffs, for cached data
    mutable std::vector<int> shCoeffIndex;
    mutable std::vector<int> shCoeffLength;
    mutable std::vector<double> shCoeffValues;

    mutable std::vector<NearEvalSH> shSrc;
    mutable std::vector<NearEvalSH> shTrg;
    mutable std::shared_ptr<InteractionManager<double, 3, NearEvalSH, NearEvalSH>> interactManagerPtr; // near neighbor
    mutable std::shared_ptr<NearInteraction<double, 3>> nearInteractorPtr;

    // dof data
    // DofIndex: beginning index of each SphHarm
    // DofLength: Length of the array
    // std::vector<int> gridValueDofIndex; // dim 3
    // std::vector<int> gridValueDofLength;
    // std::vector<double> gridValues; // excluding the north and south pole

    std::vector<int> gridNumberIndex; // dim 1
    std::vector<int> gridNumberLength;
    std::vector<double> gridWeights;        // excluding the north and south pole
    std::vector<double> gridNorms;          // excluding the north and south pole
    std::vector<double> gridCoords;         // excluding the north and south pole
    std::vector<double> gridCoordsRelative; // excluding the north and south pole, relative to the sph centers

    // fmm data
    mutable std::vector<double> srcSLCoord; // should be equal to gridPoints
    mutable std::vector<double> srcDLCoord; // should be equal to gridPoints
    mutable std::vector<double> trgCoord;
    mutable std::vector<double> srcSLValue;
    mutable std::vector<double> srcDLValue;
    mutable std::vector<double> trgValue;

    // tools
    std::shared_ptr<STKFMM> fmmPtr;
    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<TMAP> sphereMapRcp;
    Teuchos::RCP<TMAP> gridValueDofMapRcp;
    Teuchos::RCP<TMAP> gridNumberMapRcp;

    // right side of linear equation
    Teuchos::RCP<TMV> rightSideRcp;

    void setupDOF();

    void setupFMM();

    void setupNearEval();

    void applyNearEval(const double cSLex, const double cTracex) const;

    void cacheSHCoeff(double *gridValues) const;

  public:
    // Constructor
    SphereSTKSHOperator(const std::vector<Sphere> *const spherePtr, const std::string &name_,
                        std::shared_ptr<STKFMM> &fmmPtr_, const double cIdentity_ = 0, const double cSL_ = 0,
                        const double cTrac_ = 0, const double cLOP_ = 0);

    // default Destructor and copy
    ~SphereSTKSHOperator() = default;

    Teuchos::RCP<const TMAP> getDomainMap() const { return gridValueDofMapRcp; }

    Teuchos::RCP<const TMAP> getRangeMap() const { return gridValueDofMapRcp; }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    void projectNullSpace(double *inPtr) const;
    // in-place operation to remove the space in grid values which are not representable by spectral coefficients
    // results are ADDED TO the values already in outPtr

    void applyP2POP(const double *inPtr, double *outPtr, double cId_, double cSL_, double cTrac_) const;
    // apply the FMM operator cId, cSL, cDL, cTrac to inPtr and
    // results are ADDED TO the values already in outPtr

    void applyLOP(const double *inPtr, double *outPtr, double cLOP_) const;
    // apply the L operator cLOP to inPtr and
    // results are ADDED TO the values already in outPtr

    template <class Fntr>
    void setupRightSide(Fntr &fntr);

    Teuchos::RCP<TMV> &getRightSide() { return rightSideRcp; }

    // get read-only reference
    const std::vector<double> &getGridWeights() const { return gridWeights; }
    const std::vector<double> &getGridCoords() const { return gridCoords; }
    const std::vector<double> &getGridNorms() const { return gridNorms; }
    const std::vector<int> &getGridNumberIndex() const { return gridNumberIndex; }
    const std::vector<int> &getGridNumberLength() const { return gridNumberLength; }
    const std::vector<Shexp> &getSH() const { return sh; }
};

#endif