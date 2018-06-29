#ifndef SPHERESTKMOBMAT_HPP_
#define SPHERESTKMOBMAT_HPP_

#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "SphereSTKSHOperator.hpp"

#include "STKFMM/STKFMM.h"
#include "Sphere/Sphere.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/EigenDef.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"

using namespace stkfmm;

class SphereSTKMobMat : public TOP {
  private:
    std::vector<Sphere> *const spherePtr; // read only
    const std::string name;               // which sphere harmonics to work on
    const double viscosity;

    std::shared_ptr<STKFMM> fmmPtr;

    Teuchos::RCP<const TCOMM> commRcp;
    Teuchos::RCP<TMAP> sphereMapRcp;
    Teuchos::RCP<TMAP> mobMapRcp;

    // linear problem Ax = b
    Teuchos::RCP<SphereSTKSHOperator> AOpRcp;
    Teuchos::RCP<TMV> xRcp;
    Teuchos::RCP<TMV> bRcp;
    Teuchos::RCP<TMV> xLastRcp; // solution of last time

    Teuchos::RCP<Belos::SolverManager<TOP::scalar_type, TMV, TOP>> solverRcp;
    Teuchos::RCP<Belos::LinearProblem<::TOP::scalar_type, TMV, TOP>> problemRcp;

    // temporary data
    mutable std::vector<double> force;
    mutable std::vector<double> vel;
    mutable std::vector<double> rho;
    mutable std::vector<double> b;

    void testOperator();

  public:
    // Constructor
    SphereSTKMobMat(std::vector<Sphere> *const spherePtr, const std::string name_, std::shared_ptr<STKFMM> &fmmPtr_,
                    const double viscosity_);

    // default Destructor and copy
    ~SphereSTKMobMat() = default;

    Teuchos::RCP<const TMAP> getDomainMap() const { return mobMapRcp; }

    Teuchos::RCP<const TMAP> getRangeMap() const { return mobMapRcp; }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    void solveMob(const double *forcePtr, double *velPtr) const;

    const std::vector<double> &getDensitySolution() const { return rho; }

    void writeBackDensitySolution();
};

#endif