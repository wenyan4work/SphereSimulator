#ifndef SPHOPERATOR_HPP
#define SPHOPERATOR_HPP

#include "SPHSphereSystem.hpp"
#include "TpetraUtil.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

class SPHSphereOperator : public TOP {
    // 
  private:
    SPHSphereSystem *systemPtr;
    RCP<const TMAP> dofMapRcp;

  public:
    // Constructor
    explicit SPHSphereOperator(SPHSphereSystem *systemPtr_) : systemPtr{systemPtr_} {
        // tasks: 
        // 1 setup the operator
        // 2 setup the right side
        // 3 setup the iterative solver
        // 4 (optional) setup the preconditioner for the operator
        
    }

    // Destructor
    ~SPHSphereOperator() = default;

    // forbid copy
    SPHSphereOperator(const SPHSphereOperator &) = delete;
    SPHSphereOperator &operator=(const SPHSphereOperator &) = delete;

    RCP<const TMAP> getDomainMap() const {
        return dofMapRcp;
    }

    RCP<const TMAP> getRangeMap() const {
        return dofMapRcp;
    }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    // EQ 22 in note.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        // Y := beta*Y + alpha*Op(A)*X
        assert(mode == Teuchos::NO_TRANS);

        printf("Operator Applied\n");
    }
};

class LinearCombinationOperator : public TOP {
  private:
    // op1 and op2 must have the same domain and range map
    const RCP<const TOP> op1Rcp;
    const RCP<const TOP> op2Rcp;
    const double alpha, beta, gamma;

  public:
    // Constructor
    // combine alpha*op1 + beta*op2 + gamma*identity
    LinearCombinationOperator(const RCP<const TOP> &op1Rcp_, const RCP<const TOP> &op2Rcp_, double alpha_, double beta_,
                              double gamma_)
        : op1Rcp(op1Rcp_), op2Rcp(op2Rcp_), alpha(alpha_), beta(beta_), gamma(gamma_) {
        assert(op1Rcp->getDomainMap() == op2Rcp->getDomainMap());
        assert(op1Rcp->getRangeMap() == op2Rcp->getRangeMap());
    }

    // Destructor
    ~LinearCombinationOperator() = default;

    // forbid copy
    LinearCombinationOperator(const LinearCombinationOperator &) = delete;
    LinearCombinationOperator &operator=(const LinearCombinationOperator &) = delete;

    Teuchos::RCP<const TMAP> getDomainMap() const { return op1Rcp->getDomainMap(); }

    Teuchos::RCP<const TMAP> getRangeMap() const { return op1Rcp->getRangeMap(); }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    // EQ 22 in note.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        // Y := beta*Y + alpha*Op(A)*X
        assert(mode == Teuchos::NO_TRANS);
        const int nCol = X.getNumVectors();
        assert(nCol == Y.getNumVectors());

        // TODO: using Kokkos interface

        printf("Operator Applied\n");
    }
};

#endif