#include "TpetraUtil.hpp"
#include "SPHSphereSystem.hpp"

class HydroSphereOperator : public TOP {

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

class MobilityOperator : public TOP {

  private:
    // This is an implementation detail; users don't need to see it.
    HydroSphereSystem *systemPtr;
    const double Tscale;   // time scale
    const double Lscale;   // length scale
    const double Escale;   // Energy scale
    const double Fscale;   // force scale
    const double Torscale; // torque scale

  public:
    // Constructor
    MobilityOperator(HydroSphereSystem *systemPtr_, const double Tscale_, const double Lscale_, const double Escale_)
        : systemPtr{ systemPtr_ }, Tscale{ Tscale_ }, Lscale{ Lscale_ }, Escale{ Escale_ }, Fscale{ Escale_ / Lscale_ },
          Torscale{ Escale_ } {
        systemPtr->getReadyForSolve(HydroSphereSystem::SolverType::IMPLICIT);
    };

    // Destructor
    ~MobilityOperator() = default;

    // forbid copy
    MobilityOperator(const MobilityOperator &) = delete;
    MobilityOperator &operator=(const MobilityOperator) = delete;

    Teuchos::RCP<const TMAP> getDomainMap() const {
        // Get the domain Map of this Operator subclass.
        return systemPtr->getMobilityMap();
    }
    Teuchos::RCP<const TMAP> getRangeMap() const {
        // Get the range Map of this Operator subclass.
        return systemPtr->getMobilityMap();
    }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    // EQ 22 in note.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        // Y := beta*Y + alpha*Op(A)*X
        assert(mode == Teuchos::NO_TRANS);
        assert(X.getNumVectors() == Y.getNumVectors());
        assert(X.getMap()->isSameAs(*(Y.getMap())));

        auto x_2d = X.getLocalView<Kokkos::HostSpace>();
        auto y_2d = Y.getLocalView<Kokkos::HostSpace>();
        Y.modify<Kokkos::HostSpace>();
        auto &sphereIO = systemPtr->sphereIO;

        const double FscaleInv = 1 / Fscale;
        const double TorscaleInv = 1 / Torscale;

        for (int c = 0; c < x_2d.dimension_1(); c++) {
            const int sphereNumber = sphereIO.size();
// step 1, copy x to sphereIO with unit scaling
#pragma omp parallel for
            for (int i = 0; i < sphereNumber; i++) {
                sphereIO[i].force[0] = x_2d(6 * i, c) * FscaleInv;
                sphereIO[i].force[1] = x_2d(6 * i + 1, c) * FscaleInv;
                sphereIO[i].force[2] = x_2d(6 * i + 2, c) * FscaleInv;
                sphereIO[i].torque[0] = x_2d(6 * i + 3, c) * TorscaleInv;
                sphereIO[i].torque[1] = x_2d(6 * i + 4, c) * TorscaleInv;
                sphereIO[i].torque[2] = x_2d(6 * i + 5, c) * TorscaleInv;
            }

            // step 2, solve implicitly
            // Done: split the solveFullImplicit() into prepare() and solve(),
            // put prepare() into the constructor
            systemPtr->solveFullImplicit();

            // step 3, copy sphereIO to y
            const double Uscale = Lscale / Tscale;
            const double Omegascale = 1 / Tscale;
#pragma omp parallel for
            for (int i = 0; i < sphereNumber; i++) {
                const Evec3 omega = sphereIO[i].direction.cross(sphereIO[i].tdot);
                y_2d(6 * i, c) = beta * y_2d(6 * i, c) + alpha * sphereIO[i].xdot[0] * Uscale;
                y_2d(6 * i + 1, c) = beta * y_2d(6 * i + 1, c) + alpha * sphereIO[i].xdot[1] * Uscale;
                y_2d(6 * i + 2, c) = beta * y_2d(6 * i + 2, c) + alpha * sphereIO[i].xdot[2] * Uscale;
                y_2d(6 * i + 3, c) = beta * y_2d(6 * i + 3, c) + alpha * omega[0] * Omegascale;
                y_2d(6 * i + 4, c) = beta * y_2d(6 * i + 4, c) + alpha * omega[1] * Omegascale;
                y_2d(6 * i + 5, c) = beta * y_2d(6 * i + 5, c) + alpha * omega[2] * Omegascale;
            }
        }

        std::cout << "Operator applied" << std::endl;
    }
};