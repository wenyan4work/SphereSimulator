#ifndef CPSOLVER_HPP_
#define CPSOLVER_HPP_

#include "Trilinos/TpetraUtil.hpp"

#include <array>
#include <deque>
#include <vector>

using IteHistory = std::deque<std::array<double, 6>>;

class CPMatOp : public TOP {
  public:
    CPMatOp(Teuchos::RCP<TOP> mobRcp_, Teuchos::RCP<TOP> fcTransRcp_) : mobRcp(mobRcp_), fcTransRcp(fcTransRcp_) {
        // check map
        TEUCHOS_TEST_FOR_EXCEPTION(!(mobRcp->getRangeMap()->isSameAs(*(fcTransRcp->getDomainMap()))),
                                   std::invalid_argument, "Mob and Fc Maps not compatible.");
        this->forceVecRcp = Teuchos::rcp(new TV(mobRcp->getRangeMap().getConst(), true));
        this->velVecRcp = Teuchos::rcp(new TV(mobRcp->getRangeMap().getConst(), true));
    }

    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
                                   "X and Y do not have the same numbers of vectors (columns).");
        TEUCHOS_TEST_FOR_EXCEPTION(!X.getMap()->isSameAs(*Y.getMap()), std::invalid_argument,
                                   "X and Y do not have the same Map.");
        const int numVecs = X.getNumVectors();
        for (int i = 0; i < numVecs; i++) {
            const Teuchos::RCP<const TV> XcolRcp = X.getVector(i);
            Teuchos::RCP<TV> YcolRcp = Y.getVectorNonConst(i);
            // step 1 force=Fc * Xcol
            fcTransRcp->apply(*XcolRcp, *forceVecRcp, Teuchos::TRANS);
            // step 2 vel = mob * Force
            mobRcp->apply(*forceVecRcp, *velVecRcp);
            // step 3 Ycol = Fc^T * vel
            fcTransRcp->apply(*velVecRcp, *YcolRcp);
        }
    }

    Teuchos::RCP<const TMAP> getDomainMap() const {
        return this->fcTransRcp->getRangeMap(); // Get the domain Map of this Operator subclass.
    }
    Teuchos::RCP<const TMAP> getRangeMap() const {
        return this->fcTransRcp->getRangeMap(); // Get the range Map of this Operator subclass.
    }

    bool hasTransposeApply() const { return false; }

    Teuchos::RCP<TOP> mobRcp;
    Teuchos::RCP<TOP> fcTransRcp;
    Teuchos::RCP<TV> forceVecRcp;
    Teuchos::RCP<TV> velVecRcp;
};

class CPSolver {
  public:
    CPSolver(const Teuchos::RCP<const TOP> &, const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TMAP> &,
             const Teuchos::RCP<const TCOMM> &);
    CPSolver(int);

    // Nesterov Acceleration
    int LCP_APGD(Teuchos::RCP<TV> &, const double, const int, IteHistory &) const;

    // Barzilai-Borwein step length
    int LCP_BBPGD(Teuchos::RCP<TV> &, const double, const int, IteHistory &) const;

    // Minimum-Map Newtom
    int LCP_mmNewton(Teuchos::RCP<TV> &, const double, const int, IteHistory &) const;

    // Test driver
    int test_LCP(double, int, int);

  private:
    Teuchos::RCP<const TOP> ARcp;
    Teuchos::RCP<const TV> bRcp;
    Teuchos::RCP<const TMAP> mapRcp;
    Teuchos::RCP<const TCOMM> commRcp;
    int localSize;
    int globalSize;
    int myRank;
    int numProcs;
    int globalIndexMinOnLocal;
    int globalIndexMaxOnLocal;

    // functions for internal use
    void clipZero(Teuchos::RCP<TV> &) const;
    void maxXY(const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TV> &, const Teuchos::RCP<TV> &) const;
    void minXY(const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TV> &, const Teuchos::RCP<TV> &) const;
    void hMinMap(const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TV> &, const Teuchos::RCP<TV> &hRcp,
                 const Teuchos::RCP<TV> &) const;
    double checkResiduePhi(const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TV> &,
                           const Teuchos::RCP<const TV> &, const Teuchos::RCP<TV> &) const;
    double checkResiduePhi(const Teuchos::RCP<const TV> &, const Teuchos::RCP<const TV> &,
                           const Teuchos::RCP<TV> &) const;
};

#endif