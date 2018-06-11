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

class SphereSPHOperator : public TOP {
  private:
    const double cSL;
    const double cDL;
    const double cTrac; // this is ignored if the sph is a laplace layer
    const int dimension;

    const std::vector<Sphere> *const spherePtr; // read only
    const std::string name;                     // which sphere harmonics to work on

    // temporary data
    std::vector<Shexp> sph;
    // std::vector<int> spectralDofIndex;
    // std::vector<int> spectralDofOffset;

    // dof data
    std::vector<int> gridValueDofIndex;
    std::vector<int> gridValueDofOffset;
    std::vector<double> gridValues; // excluding the north and south pole

    std::vector<int> gridWeightDofIndex;
    std::vector<int> gridWeightDofOffset;
    std::vector<double> gridWeights; // excluding the north and south pole
    std::vector<double> gridPoints;  // excluding the north and south pole

    // fmm data
    std::vector<double> srcSLCoord; // should be equal to gridPoints
    std::vector<double> srcDLCoord; // should be equal to gridPoints
    std::vector<double> trgCoord;
    std::vector<double> srcSLValue;
    std::vector<double> srcDLValue;
    std::vector<double> trgValue;

    // tools
    std::shared_ptr<STKFMM> fmmPtr;
    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<TMAP> sphereMapRcp;
    Teuchos::RCP<TMAP> spectralDofMapRcp;
    Teuchos::RCP<TMAP> gridValueDofMapRcp;
    Teuchos::RCP<TMAP> gridWeightDofMapRcp;

    // right side of linear equation
    Teuchos::RCP<TMV> rightSideRcp;

    void setupDOF();

    void setupFMM();

    template <class Fntr>
    void setupRightSide(Fntr &fntr);

    Teuchos::RCP<TMV> &getRightSide() { return rightSideRcp; }

  public:
    // Constructor
    SphereSPHOperator(const std::vector<Sphere> &sphere, const std::string &name_, std::shared_ptr<STKFMM> &fmmPtr_,
                      const double cIdentity_, const double cSL_ = 0, const double cDL_ = 0, const double cTrac_ = 0);

    // Destructor
    ~SphereSPHOperator() = default;

    // forbid copy
    SphereSPHOperator(const SphereSPHOperator &) = delete;
    SphereSPHOperator &operator=(const SphereSPHOperator &) = delete;

    Teuchos::RCP<const TMAP> getDomainMap() const { return spectralDofMapRcp; }

    Teuchos::RCP<const TMAP> getRangeMap() const { return spectralDofMapRcp; }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {}
};

#endif