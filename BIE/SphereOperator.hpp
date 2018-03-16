#ifndef SPHEREBIEOPERATOR_HPP
#define SPHEREBIEOPERATOR_HPP


#include "Trilinos/TpetraUtil.hpp"
#include "Sphere/Sphere.hpp"
#include "STKFMM/STKFMM.h"

class SphereBIEOperator : public TOP {
  private:
    std::vector<Sphere> * spherePtr;
    Teuchos::RCP<TMAP> sphereMap;
    Teuchos::RCP<TMAP> dofMap;

  public:
    // Constructor
    explicit SphereBIEOperator(SPHSphereSystem *systemPtr_) : systemPtr{systemPtr_} {
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


#endif