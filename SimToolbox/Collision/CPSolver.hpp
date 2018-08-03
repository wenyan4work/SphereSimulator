#ifndef CPSOLVER_HPP_
#define CPSOLVER_HPP_

#include <array>
#include <deque>
#include <vector>

#include "Trilinos/TpetraUtil.hpp"

using IteHistory = std::deque<std::array<double, 6>>;

class CPMatOp : public TOP {
  public:
    CPMatOp(Teuchos::RCP<TOP> mobRcp_, Teuchos::RCP<TCMAT> fcTransRcp_) : mobRcp(mobRcp_), fcTransRcp(fcTransRcp_) {
        // check map
        TEUCHOS_TEST_FOR_EXCEPTION(!(mobRcp->getRangeMap()->isSameAs(*(fcTransRcp->getDomainMap()))),
                                   std::invalid_argument, "Mob and Fc Maps not compatible.");
        this->forceVecRcp = Teuchos::rcp(new TV(mobRcp->getRangeMap().getConst(), true));
        this->velVecRcp = Teuchos::rcp(new TV(mobRcp->getRangeMap().getConst(), true));
        // explicit transpose
        // Tpetra::RowMatrixTransposer<TCMAT::scalar_type, TCMAT::local_ordinal_type, TCMAT::global_ordinal_type>
        // transposer(fcTransRcp); fcRcp = transposer.createTranspose();
    }

    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
#ifdef DEBUGLCPCOL
        TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
                                   "X and Y do not have the same numbers of vectors (columns).");
        TEUCHOS_TEST_FOR_EXCEPTION(!X.getMap()->isSameAs(*Y.getMap()), std::invalid_argument,
                                   "X and Y do not have the same Map.");
#endif
        const int numVecs = X.getNumVectors();
        for (int i = 0; i < numVecs; i++) {
            const Teuchos::RCP<const TV> XcolRcp = X.getVector(i);
            Teuchos::RCP<TV> YcolRcp = Y.getVectorNonConst(i);
            // step 1 force=Fc * Xcol
            // Teuchos::RCP<Teuchos::Time> transTimer = Teuchos::TimeMonitor::getNewCounter("BBPGD::OP::FcTrans Apply");
            // {
            //     Teuchos::TimeMonitor mon(*transTimer);
            //     fcTransRcp->apply(*XcolRcp, *forceVecRcp, Teuchos::TRANS);
            //     // fcRcp->apply(*XcolRcp, *forceVecRcp);
            // }
            // // step 2 vel = mob * Force
            // Teuchos::RCP<Teuchos::Time> mobTimer = Teuchos::TimeMonitor::getNewCounter("BBPGD::OP::Mob Apply");
            // {
            //     Teuchos::TimeMonitor mon(*mobTimer);
            //     mobRcp->apply(*forceVecRcp, *velVecRcp);
            // }
            // // step 3 Ycol = Fc^T * vel
            // Teuchos::RCP<Teuchos::Time> fcTimer = Teuchos::TimeMonitor::getNewCounter("BBPGD::OP::Fc Apply");
            // {
            //     Teuchos::TimeMonitor mon(*fcTimer);
            //     fcTransRcp->apply(*velVecRcp, *YcolRcp);
            // }
            fcTransRcp->apply(*XcolRcp, *forceVecRcp, Teuchos::TRANS);
            mobRcp->apply(*forceVecRcp, *velVecRcp);
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
    Teuchos::RCP<TCMAT> fcTransRcp;
    // Teuchos::RCP<TCMAT> fcRcp;
    Teuchos::RCP<TV> forceVecRcp;
    Teuchos::RCP<TV> velVecRcp;
};

class CPSolver {
  public:
    CPSolver(const Teuchos::RCP<const TOP> &A_, const Teuchos::RCP<const TV> &b_);
    CPSolver(int localSize, double diagonal = 0.0);

    // Nesterov Acceleration
    int LCP_APGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const;

    // Barzilai-Borwein step length
    int LCP_BBPGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const;

    // Minimum-Map Newtom
    int LCP_mmNewton(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const;

    // Test driver
    int test_LCP(double tol, int maxIte, int solverChoice);

  private:
    Teuchos::RCP<const TOP> ARcp;
    Teuchos::RCP<const TV> bRcp;
    Teuchos::RCP<const TMAP> mapRcp; // map for the distribution of xsolRcp / bRcp / ARcp->rowMap
    Teuchos::RCP<const TCOMM> commRcp;

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