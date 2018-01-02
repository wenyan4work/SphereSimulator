#include "TpetraUtil.hpp"

class SphereOperator : public TOP {

  private:
    // This is an implementation detail; users don't need to see it.
    HydroSphereSystem *systemPtr;

  public:
    // Constructor
    explicit HydroSphereOperator(HydroSphereSystem *systemPtr_) : systemPtr{ systemPtr_ } {}

    // Destructor
    ~HydroSphereOperator() = default;

    // forbid copy
    HydroSphereOperator(const HydroSphereOperator &) = delete;
    HydroSphereOperator &operator=(const HydroSphereOperator) = delete;

    Teuchos::RCP<const TMAP> getDomainMap() const {
        // Get the domain Map of this Operator subclass.
        return systemPtr->getFdistMap();
    }

    Teuchos::RCP<const TMAP> getRangeMap() const {
        // Get the range Map of this Operator subclass.
        return systemPtr->getFdistMap();
    }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    // EQ 22 in note.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        // Y := beta*Y + alpha*Op(A)*X
        assert(mode == Teuchos::NO_TRANS);
        systemPtr->applyImpSolverMatrixFdist(X, Y, alpha, beta);
        std::cout << "Operator applied" << std::endl;
    }
};
