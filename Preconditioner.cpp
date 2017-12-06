#include "Preconditioner.hpp"

void showPrecDebugInfo(const Teuchos::RCP<const TCMAT> &A, const Teuchos::RCP<const TOP> &PrecOp) {
#ifdef IFPACKDEBUG
    dumpTCMAT(A, "A_for_Prec.mtx");
    if (A->getComm()->getRank() == 0) {
        Teuchos::FancyOStream fancyout(Teuchos::rcp(&std::cout, false));
        cout << "NNZ local" << A->getNodeNumEntries() << endl;

        cout << "A domainMap" << A->getDomainMap()->description() << endl;
        cout << "A domainMap Min Local" << A->getDomainMap()->getMinLocalIndex() << endl;
        cout << "A domainMap Max Local" << A->getDomainMap()->getMaxLocalIndex() << endl;
        cout << "A domainMap Min Global" << A->getDomainMap()->getMinGlobalIndex() << endl;
        cout << "A domainMap Max Global" << A->getDomainMap()->getMaxGlobalIndex() << endl;

        cout << "A rangeMap" << A->getRangeMap()->description() << endl;
        cout << "A rangeMap Min Local" << A->getRangeMap()->getMinLocalIndex() << endl;
        cout << "A rangeMap Max Local" << A->getRangeMap()->getMaxLocalIndex() << endl;
        cout << "A rangeMap Min Global" << A->getRangeMap()->getMinGlobalIndex() << endl;
        cout << "A rangeMap Max Global" << A->getRangeMap()->getMaxGlobalIndex() << endl;

        cout << "A rowMap" << A->getRowMap()->description() << endl;
        cout << "A rowMap Min Local" << A->getRowMap()->getMinLocalIndex() << endl;
        cout << "A rowMap Max Local" << A->getRowMap()->getMaxLocalIndex() << endl;
        cout << "A rowMap Min Global" << A->getRowMap()->getMinGlobalIndex() << endl;
        cout << "A rowMap Max Global" << A->getRowMap()->getMaxGlobalIndex() << endl;

        cout << "A colMap" << A->getColMap()->description() << endl;
        cout << "A colMap Min Local" << A->getColMap()->getMinLocalIndex() << endl;
        cout << "A colMap Max Local" << A->getColMap()->getMaxLocalIndex() << endl;
        cout << "A colMap Min Global" << A->getColMap()->getMinGlobalIndex() << endl;
        cout << "A colMap Max Global" << A->getColMap()->getMaxGlobalIndex() << endl;
    }
    A->getComm()->barrier();
    if (A->getComm()->getRank() == 0) {
        cout << "description" << prec->description() << endl;
        Teuchos::FancyOStream fancyout(Teuchos::rcp(&std::cout, false));
        prec->describe(fancyout, Teuchos::EVerbosityLevel::VERB_EXTREME);
        prec->getDomainMap()->describe(fancyout, Teuchos::EVerbosityLevel::VERB_EXTREME);
        prec->getRangeMap()->describe(fancyout, Teuchos::EVerbosityLevel::VERB_EXTREME);
        cout << "prec domainMap" << prec->getDomainMap()->description() << endl;
        cout << "prec domainMap Min Local" << prec->getDomainMap()->getMinLocalIndex() << endl;
        cout << "prec domainMap Max Local" << prec->getDomainMap()->getMaxLocalIndex() << endl;
        cout << "prec domainMap Min Global" << prec->getDomainMap()->getMinGlobalIndex() << endl;
        cout << "prec domainMap Max Global" << prec->getDomainMap()->getMaxGlobalIndex() << endl;
    }
    Teuchos::RCP<TMV> xRCP = Teuchos::rcp(new TMV(A->getRowMap(), 1, false));
    Teuchos::RCP<TMV> bRCP = Teuchos::rcp(new TMV(A->getRowMap(), 1, false));
    xRCP->putScalar(1);
    A->apply(*xRCP, *bRCP);
    dumpTMV(bRCP, "bRCP_A.mtx");
    prec->apply(*xRCP, *bRCP);
    dumpTMV(bRCP, "bRCP_prec,mtx");
#endif
}

// ILUT preconditioner, not good for threading, and includes only diagonal blocks in Schwarz decomposition
Teuchos::RCP<TOP> createILUTPreconditioner(const Teuchos::RCP<const TCMAT> &A, double tol, double fill) {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using std::cout;
    using std::endl;

    // set parameter list
    // The name of the type of preconditioner to use.
    Teuchos::ParameterList plist;

    // parameters for ILUT
    const std::string precondType = "ILUT";
    plist.set("fact: ilut level-of-fill", fill);
    plist.set("fact: drop tolerance", tol);
    plist.set("fact: absolute threshold", 0.0001); // this stablizes the preconditioner but takes more iterations

    // Fetch the typedefs defined by Tpetra::CrsMatrix.
    typedef typename TCMAT::scalar_type scalar_type;
    typedef typename TCMAT::local_ordinal_type local_ordinal_type;
    typedef typename TCMAT::global_ordinal_type global_ordinal_type;
    //	typedef typename TpetraMatrixType::node_type node_type;

    typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type> op_type;

    // These are just some convenience typedefs.
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    // An Ifpack2::Preconditioner is-a Tpetra::Operator.  Ifpack2
    // creates a Preconditioner object, but users of iterative methods
    // want a Tpetra::Operator.  That's why create() returns an Operator
    // instead of a Preconditioner.
    typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type> prec_type;

    // Create timers to show how long it takes for Ifpack2 to do various operations.
    RCP<Time> initTimer = TimeMonitor::getNewCounter("Ifpack2::Preconditioner::initialize");
    RCP<Time> computeTimer = TimeMonitor::getNewCounter("Ifpack2::Preconditioner::compute");
    RCP<Time> condestTimer = TimeMonitor::getNewCounter("Ifpack2::Preconditioner::condest");

    cout << "Creating ILUT preconditioner\n "
         << "-- Configuring" << endl;
    //
    // Create the preconditioner and set parameters.
    //
    // This doesn't actually _compute_ the preconditioner.
    // It just sets up the specific type of preconditioner and
    // its associated parameters (which depend on the type).
    RCP<prec_type> prec;
    Ifpack2::Factory factory;
    // Set up the preconditioner of the given type.
    prec = factory.create(precondType, A);
    prec->setParameters(plist);

    cout << "-- Initializing" << endl;
    {
        TimeMonitor mon(*initTimer);
        prec->initialize();
    }

    // THIS ACTUALLY COMPUTES THE PRECONDITIONER
    // (e.g., does the incomplete factorization).
    cout << "-- Computing" << endl;
    {
        TimeMonitor mon(*computeTimer);
        prec->compute();
    }

    showPrecDebugInfo(A, prec);

    return prec;
}

// Polynomial Relaxation Preconditioner
Teuchos::RCP<TOP> createPlnPreconditioner(const Teuchos::RCP<const TCMAT> &A) {
    // get a CrsMatrix A, return a Polynomial (Relaxation, Chebyshev, etc) Preconditioner for A represented by a
    // Operator
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using std::endl;
    using std::cout;
    // set parameter list
    // The name of the type of preconditioner to use.
    //	Teuchos::ParameterList plistSWZ; // SCHWARZ domain decomposition -> block diagonal
    //	const std::string precondType = "SCHWARZ";
    const std::string precondType = "RELAXATION";
    //	plistSWZ.set("schwarz: inner preconditioner name", "ILUT");
    //	plistSWZ.set("schwarz: overlap level", 50); set this to non zero greatly deteriorates the iterations.
    //	plistSWZ.set("schwarz: combine mode", "ADD");

    Teuchos::ParameterList plist; // method for each local diagonal block
                                  //	plist.set("fact: ilut level-of-fill",5);
                                  //	plist.set("fact: drop tolerance", 1e-3);
                                  //	plist.set("fact: absolute threshold", 1e-5);

    plist.set("relaxation: type", "Jacobi"); // relaxation effective only for intermediate volume fraction without
                                             // Schwarz. Jacobi shows good threading in this case
    plist.set("relaxation: sweeps", 5);      // further increase sweeps is probably useless
    plist.set("relaxation: damping factor", 1.0);
    // reference: for two fibers close to each other
    // Jacobi sweep = 5:
    // damping 1 -> 2 iterations
    // damping 2 -> 9 iterations
    // damping 0.5 -> 4  iterations
    // sweep 10, damping 0.8 -> 2 ite
    // sweep 5, damping 0.8 -> 4 ite
    // Gauss-Seidel sweep = 5:
    // damping 0.5 -> 8 iterations
    // damping 1.0 -> 9 iterations
    // damping 2.0 -> 40 iterations
    // Symmetric Gauss-Seidel sweep = 5:
    // damping 0.5 -> 8 iterations
    // damping 1.0 -> 19 iterations
    // damping 2.0 -> too many, cannot finish iterations
    plist.set("relaxation: use l1", true); // a fix for MPI?
    plist.set("relaxation: fix tiny diagonal entries", true);
    plist.set("relaxation: zero starting solution", true); // must be true. otherwise may give random or NAN result
    plist.set("relaxation: min diagonal value", 1e-5);

    //	plistSWZ.set("schwarz: inner preconditioner parameters", plist);

    // Fetch the typedefs defined by Tpetra::CrsMatrix.
    typedef typename TCMAT::scalar_type scalar_type;
    typedef typename TCMAT::local_ordinal_type local_ordinal_type;
    typedef typename TCMAT::global_ordinal_type global_ordinal_type;
    //	typedef typename TpetraMatrixType::node_type node_type;

    // Ifpack2's generic Preconditioner interface implements
    // Tpetra::Operator.  A Tpetra::Operator is an abstraction of a
    // function mapping a (Multi)Vector to a (Multi)Vector, with the
    // option of applying the transpose or conjugate transpose of the
    // operator.  Tpetra::CrsMatrix implements Operator as well.
    typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type> op_type;

    // These are just some convenience typedefs.
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef typename STS::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;

    typedef Ifpack2::Preconditioner<scalar_type, local_ordinal_type, global_ordinal_type> prec_type;

    // Create timers to show how long it takes for Ifpack2 to do various operations.
    RCP<Time> initTimer = TimeMonitor::getNewCounter("Ifpack2::Preconditioner::initialize");
    RCP<Time> computeTimer = TimeMonitor::getNewCounter("Ifpack2::Preconditioner::compute");

    cout << "Creating preconditioner\n"
         << "-- Configuring" << endl;

    // Create the preconditioner and set parameters.
    //
    // This doesn't actually _compute_ the preconditioner.
    // It just sets up the specific type of preconditioner and
    // its associated parameters (which depend on the type).
    RCP<prec_type> prec;
    Ifpack2::Factory factory;
    // Set up the preconditioner of the given type.
    prec = factory.create(precondType, A);
    prec->setParameters(plist);
    cout << "-- Initializing" << endl;
    {
        TimeMonitor mon(*initTimer);
        prec->initialize();
    }

    // THIS ACTUALLY COMPUTES THE PRECONDITIONER
    // (e.g., does the incomplete factorization).
    cout << "-- Computing" << endl;
    {
        TimeMonitor mon(*computeTimer);
        prec->compute();
    }

    showPrecDebugInfo(A, prec);
    return prec;
}

class KinvOperator : public TOP {
  private:
    // This is an implementation detail; users don't need to see it.
    const Teuchos::RCP<const TCMAT> A;
    const Teuchos::RCP<const Teuchos::Comm<int>> comm;
    const Teuchos::RCP<const TMAP> rowMap;
    const Teuchos::RCP<const TMAP> colMap;
    const int localSize;
    const int globalSize;
    const int myRank;
    const int numProcs;
    const int globalIndexMinOnLocal;
    const int globalIndexMaxOnLocal;
    // set Belos object
    Belos::SolverFactory<::TOP::scalar_type, TMV, TOP> factory;
    Teuchos::RCP<Belos::SolverManager<::TOP::scalar_type, TMV, TOP>> solverRCP;
    Teuchos::RCP<Belos::LinearProblem<::TOP::scalar_type, TMV, TOP>> problemRCP;
    Teuchos::RCP<TMV> inMV;
    Teuchos::RCP<TMV> outMV;
    Teuchos::RCP<TOP> precOp;

  public:
    // Constructor
    //
    // n: Global number of rows and columns in the operator.
    // comm: The communicator over which to distribute those rows and columns.
    KinvOperator(const Teuchos::RCP<const TCMAT> &A, const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                 const Teuchos::RCP<const TMAP> &map, const Teuchos::RCP<const TMV> &initialGuess)
        // initialGuess vector is not modified.
        : KinvOperator(A, initialGuess){}; // delegate constructor

    KinvOperator(const Teuchos::RCP<const TCMAT> &A_, const Teuchos::RCP<const TMV> &initialGuess_)
        : A(A_), comm(A_->getComm()), rowMap(A_->getRowMap()), colMap(A_->getColMap()), myRank(comm->getRank()),
          numProcs(comm->getSize()), globalSize(rowMap->getGlobalNumElements()),
          localSize(rowMap->getNodeNumElements()), globalIndexMinOnLocal(rowMap->getMinGlobalIndex()),
          globalIndexMaxOnLocal(rowMap->getMaxGlobalIndex()) {
        TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument,
                                   "MyOp constructor: The input Comm object must be nonnull.");

        inMV = Teuchos::rcp(new TMV(*initialGuess_, Teuchos::Copy));
        outMV = Teuchos::rcp(new TMV(rowMap, 1, true)); // true means zero initial value
                                                        //		inMV->putScalar(1.0); // a valid initial guess
                                                        //		outMV->putScalar(1.0); // a valid initial guess

        // Make an empty new parameter list.
        Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
        solverParams->set("Num Blocks", 200); // larger than this might trigger a std::bad_alloc inside Kokkos.
        solverParams->set("Maximum Iterations", 1000);
        solverParams->set("Convergence Tolerance", 1e-8);
        solverParams->set("Timer Label", "Iterative Inverse Preconditioner");
#ifdef IFPACKDEBUG
        solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);
#endif
        // Create the GMRES solver.
        solverRCP = factory.create("GMRES", solverParams); // GCRODR does not work well in this case
        //		solverRCP->reset(Belos::Problem);
        //		solverRCP->reset(Belos::RecycleSubspace);

        // get prec
        comm->barrier();
        precOp = createPlnPreconditioner(A.getConst());
        // rightPrecondition with this in GCRODR has a bug of MV dimension

        // setup the problem
        problemRCP = Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(A, outMV, inMV));
        problemRCP->setProblem();
        problemRCP->setRightPrec(precOp);
        solverRCP->setProblem(problemRCP);
        comm->barrier();
#ifdef IFPACKDEBUG
        if (myRank == 0) {
            std::cout << precOp->description() << std::endl;
            std::cout << "KinvOperator constructed" << solverRCP->description() << std::endl;
        }
#endif
    }

    // Destructor
    ~KinvOperator() = default;

    Teuchos::RCP<const TMAP> getDomainMap() const {
        return this->rowMap; // Get the domain Map of this Operator subclass.
    }
    Teuchos::RCP<const TMAP> getRangeMap() const {
        return this->rowMap; // Get the range Map of this Operator subclass.
    }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    // EQ 22 in note.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        Teuchos::RCP<Teuchos::Time> precTimer = Teuchos::TimeMonitor::getNewCounter("Preconditioner Solve");

        {
            Teuchos::TimeMonitor mon(*precTimer);
            inMV->update(1.0, X, 0.0);
            outMV->update(1.0, X, 0.0);

            solverRCP->reset(Belos::Problem);
#ifdef IFPACKDEBUG
            if (myRank == 0) {
                std::cout << solverRCP->description() << std::endl;
                solverRCP->describe(std::cout, Teuchos::VERB_EXTREME);
                std::cout << "KinvOperator reset" << std::endl << solverRCP->description() << std::endl;
            }
#endif
            Belos::ReturnType result = solverRCP->solve();
            int numIters2 = solverRCP->getNumIters();
            if (myRank == 0) {
                std::cout << "numIters Preconditioner: " << numIters2 << std::endl;
            }

            Y.update(1.0, *outMV, 0.0);
        }
    }
};

Teuchos::RCP<TOP> createKinvPreconditioner(const Teuchos::RCP<const TCMAT> &A,
                                           const Teuchos::RCP<const TMV> &initialGuess) {
    return Teuchos::rcp(new KinvOperator(A.getConst(), initialGuess.getConst()));
}
