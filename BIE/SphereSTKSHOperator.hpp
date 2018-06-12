#ifndef SPHERESPHOPERATOR_HPP
#define SPHERESPHOPERATOR_HPP

#include "STKFMM/STKFMM.h"
#include "Sphere/Sphere.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include <mpi.h>
#include <omp.h>

using namespace stkfmm;

class SphereSTKSHOperator : public TOP {
  private:
    const double cSL;
    const double cDL;
    const double cTrac;
    const double cLOP;
    const int dimension;

    const std::vector<Sphere> *const spherePtr; // read only
    const std::string name;                     // which sphere harmonics to work on

    // temporary data
    mutable std::vector<Shexp> sph;
    mutable std::vector<double> pointValues;

    // dof data
    // DofIndex: beginning index of each SphHarm
    // DofLength: Length of the array
    std::vector<int> gridValueDofIndex; // dim 3
    std::vector<int> gridValueDofLength;
    std::vector<double> gridValues; // excluding the north and south pole

    std::vector<int> gridDofIndex; // dim 1
    std::vector<int> gridDofLength;
    std::vector<double> gridWeights; // excluding the north and south pole
    std::vector<double> gridNorms;   // excluding the north and south pole
    std::vector<double> gridCoords;  // excluding the north and south pole

    // fmm data
    std::vector<double> srcSLCoord; // should be equal to gridPoints
    std::vector<double> srcDLCoord; // should be equal to gridPoints
    std::vector<double> trgCoord;
    mutable std::vector<double> srcSLValue;
    mutable std::vector<double> srcDLValue;
    mutable std::vector<double> trgValue;

    // tools
    std::shared_ptr<STKFMM> fmmPtr;
    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<TMAP> sphereMapRcp;
    Teuchos::RCP<TMAP> gridValueDofMapRcp;
    Teuchos::RCP<TMAP> gridDofMapRcp;

    // right side of linear equation
    Teuchos::RCP<TMV> rightSideRcp;

    void setupDOF();

    void setupFMM();

    template <class Fntr>
    void setupRightSide(Fntr &fntr);

    void projectNullSpace() const;

    void runFMM() const;

  public:
    // Constructor
    SphereSTKSHOperator(const std::vector<Sphere> &sphere, const std::string &name_, std::shared_ptr<STKFMM> &fmmPtr_,
                        const double cIdentity_ = 0, const double cSL_ = 0, const double cDL_ = 0,
                        const double cTrac_ = 0, const double cLOP_ = 0);

    // Destructor
    ~SphereSTKSHOperator() = default;

    // forbid copy
    SphereSTKSHOperator(const SphereSTKSHOperator &) = delete;
    SphereSTKSHOperator &operator=(const SphereSTKSHOperator &) = delete;

    Teuchos::RCP<const TMAP> getDomainMap() const { return gridValueDofMapRcp; }

    Teuchos::RCP<const TMAP> getRangeMap() const { return gridValueDofMapRcp; }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {}

    Teuchos::RCP<TMV> &getRightSide() { return rightSideRcp; }
};

#endif