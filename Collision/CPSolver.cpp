
#include <chrono>

#include "CPSolver.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "mpi.h"

void CPSolver::clipZero(Teuchos::RCP<TV> &vecRcp) const {
    auto x_2d = vecRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    vecRcp->modify<Kokkos::HostSpace>();

    // int strides[2];
    // x_2d.stride(strides);
    // std::cout << strides[0] << " " << strides[1] << std::endl;

    // std::cout << std::is_same<decltype(x_2d)::array_layout, Kokkos::LayoutLeft>::value << std::endl;

    for (int c = 0; c < x_2d.dimension_1(); c++) {
#pragma omp parallel for
        for (int i = 0; i < x_2d.dimension_0(); i++) {
            const double temp = x_2d(i, c);
            x_2d(i, c) = temp < 0 ? 0 : temp;
        }
    }
    return;
}

void CPSolver::maxXY(const Teuchos::RCP<const TV> &vecXRcp, const Teuchos::RCP<const TV> &vecYRcp,
                     const Teuchos::RCP<TV> &vecZRcp) const {
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getMap()->isSameAs(*(vecZRcp->getMap()))), std::invalid_argument,
                               "X and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getMap()->isSameAs(*(vecZRcp->getMap()))), std::invalid_argument,
                               "Y and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getNumVectors() == vecZRcp->getNumVectors()), std::invalid_argument,
                               "X and Z do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getNumVectors() == vecZRcp->getNumVectors()), std::invalid_argument,
                               "Y and Z do not have the same Number of Vectors.");
#endif

    auto x_2d = vecXRcp->getLocalView<Kokkos::HostSpace>();
    auto y_2d = vecYRcp->getLocalView<Kokkos::HostSpace>();
    auto z_2d = vecZRcp->getLocalView<Kokkos::HostSpace>();
    vecZRcp->modify<Kokkos::HostSpace>();
    for (int c = 0; c < x_2d.dimension_1(); c++) {
#pragma omp parallel for
        for (int i = 0; i < x_2d.dimension_0(); i++) {
            z_2d(i, c) = std::max(x_2d(i, c), y_2d(i, c));
        }
    }
    return;
}

void CPSolver::minXY(const Teuchos::RCP<const TV> &vecXRcp, const Teuchos::RCP<const TV> &vecYRcp,
                     const Teuchos::RCP<TV> &vecZRcp) const {
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getMap()->isSameAs(*(vecZRcp->getMap()))), std::invalid_argument,
                               "X and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getMap()->isSameAs(*(vecZRcp->getMap()))), std::invalid_argument,
                               "Y and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getNumVectors() == vecZRcp->getNumVectors()), std::invalid_argument,
                               "X and Z do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getNumVectors() == vecZRcp->getNumVectors()), std::invalid_argument,
                               "Y and Z do not have the same Number of Vectors.");
#endif

    auto x_2d = vecXRcp->getLocalView<Kokkos::HostSpace>();
    auto y_2d = vecYRcp->getLocalView<Kokkos::HostSpace>();
    auto z_2d = vecZRcp->getLocalView<Kokkos::HostSpace>();
    vecZRcp->modify<Kokkos::HostSpace>();
    for (int c = 0; c < x_2d.dimension_1(); c++) {
#pragma omp parallel for
        for (int i = 0; i < x_2d.dimension_0(); i++) {
            z_2d(i, c) = std::min(x_2d(i, c), y_2d(i, c));
        }
    }
    return;
}

double CPSolver::checkResiduePhi(const Teuchos::RCP<const TV> &vecXRcp, const Teuchos::RCP<const TV> &vecYRcp,
                                 const Teuchos::RCP<const TV> &vecbRcp, const Teuchos::RCP<TV> &vecTempRcp) const {
// Y=Ax
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getMap()->isSameAs(*(vecTempRcp->getMap()))), std::invalid_argument,
                               "X and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getMap()->isSameAs(*(vecTempRcp->getMap()))), std::invalid_argument,
                               "Y and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getNumVectors() == vecTempRcp->getNumVectors()), std::invalid_argument,
                               "X and Z do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getNumVectors() == vecTempRcp->getNumVectors()), std::invalid_argument,
                               "Y and Z do not have the same Number of Vectors.");
#endif
    auto x_2d = vecXRcp->getLocalView<Kokkos::HostSpace>();
    auto y_2d = vecYRcp->getLocalView<Kokkos::HostSpace>();
    auto b_2d = vecbRcp->getLocalView<Kokkos::HostSpace>();
    auto z_2d = vecTempRcp->getLocalView<Kokkos::HostSpace>();
    vecTempRcp->modify<Kokkos::HostSpace>();
    for (int c = 0; c < x_2d.dimension_1(); c++) {
#pragma omp parallel for
        for (int i = 0; i < x_2d.dimension_0(); i++) {
            z_2d(i, c) = std::min(x_2d(i, c), y_2d(i, c) + b_2d(i, c));
        }
    }
    return static_cast<double>(vecTempRcp->norm2());
}

double CPSolver::checkResiduePhi(const Teuchos::RCP<const TV> &vecXRcp, const Teuchos::RCP<const TV> &vecYRcp,
                                 const Teuchos::RCP<TV> &vecTempRcp) const {
// Y=Ax+b
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getMap()->isSameAs(*(vecTempRcp->getMap()))), std::invalid_argument,
                               "X and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getMap()->isSameAs(*(vecTempRcp->getMap()))), std::invalid_argument,
                               "Y and Z do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecXRcp->getNumVectors() == vecTempRcp->getNumVectors()), std::invalid_argument,
                               "X and Z do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecYRcp->getNumVectors() == vecTempRcp->getNumVectors()), std::invalid_argument,
                               "Y and Z do not have the same Number of Vectors.");
#endif
    minXY(vecXRcp, vecYRcp, vecTempRcp);
    return static_cast<double>(vecTempRcp->norm2());
}

void CPSolver::hMinMap(const Teuchos::RCP<const TV> &xRcp, const Teuchos::RCP<const TV> &yRcp,
                       const Teuchos::RCP<TV> &hRcp, const Teuchos::RCP<TV> &maskRcp) const {
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(xRcp->getMap()->isSameAs(*(hRcp->getMap()))), std::invalid_argument,
                               "X and h do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(yRcp->getMap()->isSameAs(*(hRcp->getMap()))), std::invalid_argument,
                               "Y and h do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(xRcp->getNumVectors() == hRcp->getNumVectors()), std::invalid_argument,
                               "X and h do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(yRcp->getNumVectors() == hRcp->getNumVectors()), std::invalid_argument,
                               "Y and h do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(xRcp->getMap()->isSameAs(*(maskRcp->getMap()))), std::invalid_argument,
                               "X and mask do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(yRcp->getMap()->isSameAs(*(maskRcp->getMap()))), std::invalid_argument,
                               "Y and mask do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(xRcp->getNumVectors() == maskRcp->getNumVectors()), std::invalid_argument,
                               "X and mask do not have the same Number of Vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(yRcp->getNumVectors() == maskRcp->getNumVectors()), std::invalid_argument,
                               "Y and mask do not have the same Number of Vectors.");
#endif
    // assume they are all of the same size, comm and map, do not check
    auto xView = xRcp->getLocalView<Kokkos::HostSpace>();
    auto yView = yRcp->getLocalView<Kokkos::HostSpace>();
    auto hView = hRcp->getLocalView<Kokkos::HostSpace>();
    auto maskView = maskRcp->getLocalView<Kokkos::HostSpace>();
    hRcp->modify<Kokkos::HostSpace>();
    maskRcp->modify<Kokkos::HostSpace>();

    for (int c = 0; c < xView.dimension_1(); c++) {
#pragma omp parallel for
        for (int i = 0; i < xView.dimension_0(); i++) {
            if (xView(i, c) < yView(i, c)) {
                hView(i, c) = xView(i, c);
                maskView(i, c) = 1;
            } else {
                hView(i, c) = yView(i, c);
                maskView(i, c) = 0;
            }
        }
    }

    return;
}

// constructor with known A and b
CPSolver::CPSolver(const Teuchos::RCP<const TOP> &A_, const Teuchos::RCP<const TV> &b_,
                   const Teuchos::RCP<const TMAP> &map_, const Teuchos::RCP<const TCOMM> &comm_)
    : ARcp(A_), bRcp(b_), mapRcp(map_), commRcp(comm_), myRank(comm_->getRank()), numProcs(comm_->getSize()),
      globalSize(map_->getGlobalNumElements()), localSize(map_->getNodeNumElements()),
      globalIndexMinOnLocal(map_->getMinGlobalIndex()), globalIndexMaxOnLocal(map_->getMaxGlobalIndex()) {
// make sure A and b match the map and comm specified
#ifdef DEBUGCOLLCP
    TEUCHOS_TEST_FOR_EXCEPTION(!(ARcp->getDomainMap()->isSameAs(*(bRcp->getMap()))), std::invalid_argument,
                               "A and b do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(mapRcp->isSameAs(*(bRcp->getMap()))), std::invalid_argument,
                               "map and b do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(mapRcp->isSameAs(*(ARcp->getDomainMap()))), std::invalid_argument,
                               "map and b do not have the same Map.");
#endif
}

// constructor to set random A and b with given size
CPSolver::CPSolver(int globalSize) {
    // set up comm
    Teuchos::RCP<const TCOMM> commRcp = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    myRank = commRcp->getRank();
    numProcs = commRcp->getSize();
    // set up row and col maps
    Teuchos::RCP<const TMAP> colMapRcp = Teuchos::rcp(new TMAP(globalSize, 0, commRcp));
    Teuchos::RCP<const TMAP> rowMapRcp = Teuchos::rcp(new TMAP(globalSize, 0, commRcp));
    mapRcp = rowMapRcp;
    localSize = rowMapRcp->getNodeNumElements();
    globalIndexMinOnLocal = rowMapRcp->getMinGlobalIndex();
    globalIndexMaxOnLocal = rowMapRcp->getMaxGlobalIndex();
    if (myRank == 0) {
        std::cout << "Total number of processes: " << numProcs << std::endl;
        std::cout << "rank: " << myRank << std::endl;
        std::cout << "global size: " << globalSize << std::endl;
        std::cout << "local size: " << localSize << std::endl;
        std::cout << "global index min on local: " << globalIndexMinOnLocal << std::endl;
        std::cout << "global index max on local: " << globalIndexMaxOnLocal << std::endl;
    }

    // make sure A and b match the map and comm specified
    // set A and b randomly. maintain SPD of A
    Teuchos::RCP<TV> btemp = Teuchos::rcp(new TV(rowMapRcp, false));
    btemp->randomize(-0.5, 0.5);
    bRcp = btemp.getConst();

    // generate a local random matrix
    // Create an empty matrix
    Teuchos::SerialDenseMatrix<int, double> ArootLocal(localSize, localSize, true); // zeroOut
    for (int i = 0; i < localSize; i++) {
        for (int j = 0; j < 5; j++) {
            int rowIndex = i;
            int colIndex = drand48() * localSize;
            ArootLocal(rowIndex, colIndex) = drand48() - 0.5;
        }
        ArootLocal(i, i) += 1; // add value to diagonal. Otherwise  solvers fail
    }

    Teuchos::SerialDenseMatrix<int, double> ALocal(localSize, localSize);
    ALocal.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, ArootLocal, ArootLocal, 0.0);

    // use ALocal as local matrix to fill TCMAT A
    double droptol = 1e-7;
    Kokkos::View<size_t *> rowCount("rowCount", localSize);
    Kokkos::View<size_t *> rowPointers("rowPointers", localSize + 1);
    for (int i = 0; i < localSize; i++) {
        rowCount[i] = 0;
        for (int j = 0; j < localSize; j++) {
            if (fabs(ALocal(i, j)) > droptol) {
                rowCount[i]++;
            }
        }
    }

    rowPointers[0] = 0;
    for (int i = 1; i < localSize + 1; i++) {
        rowPointers[i] = rowPointers[i - 1] + rowCount[i - 1];
    }
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localSize]);
    Kokkos::View<double *> values("values", rowPointers[localSize]);
    int p = 0;
    for (int i = 0; i < localSize; i++) {
        for (int j = 0; j < localSize; j++) {
            if (fabs(ALocal(i, j)) > droptol) {
                columnIndices[p] = j;
                values[p] = ALocal(i, j);
                p++;
            }
        }
    }

    // fill matrix Aroot
    Teuchos::RCP<TCMAT> Atemp = Teuchos::rcp(new TCMAT(rowMapRcp, colMapRcp, rowPointers, columnIndices, values));
    Atemp->fillComplete(rowMapRcp, rowMapRcp);
    this->ARcp = Teuchos::rcp_dynamic_cast<const TOP>(Atemp, true);
    std::cout << "ARcp" << ARcp->description() << std::endl;

    // dump matrix
    dumpTCMAT(Atemp, "Amat.mtx");
    dumpTV(this->bRcp, "bvec.mtx");
}

int CPSolver::LCP_BBPGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const {
    MPI_Barrier(MPI_COMM_WORLD);
    int mvCount = 0; // count matrix-vector multiplications
    if (myRank == 0) {
        std::cout << "solving BBPGD" << std::endl;
        std::cout << "ARcp" << ARcp->description() << std::endl;
    }
    // convention in the iteratin:

    // map must match
    TEUCHOS_TEST_FOR_EXCEPTION(!this->mapRcp->isSameAs(*(xsolRcp->getMap())), std::invalid_argument,
                               "xsolrcp and A operator do not have the same Map.");
    Teuchos::RCP<TV> xkm1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // zero out
    Teuchos::RCP<TV> ykm1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // zero out
    Teuchos::RCP<TV> xkRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy));         // deep copy, xk=x0
    Teuchos::RCP<TV> ykRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    Teuchos::RCP<TV> xkp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    Teuchos::RCP<TV> ykp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));

    Teuchos::RCP<TV> tempVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    Teuchos::RCP<TV> ykmykm1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    Teuchos::RCP<TV> xkmxkm1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));

    ARcp->apply(*xkm1Rcp, *ykm1Rcp); // ykm1=A.dot(xkm1)
    mvCount++;
    MPI_Barrier(MPI_COMM_WORLD);
    tempVecRcp->update(1.0, *bRcp, 1.0, *ykm1Rcp, 0.0); // tempVec=A.dot(xkm1)+b
    // first step, use the step size (x^T A g + 0.5*b^T g) / (g^T A g), global minimum of line search along g
    // use ykp1 and xkp1 as temporary space
    ARcp->apply(*tempVecRcp, *ykp1Rcp); // ykp1 =  A*tempVec = Ag
    mvCount++;
    double xTAg = xkm1Rcp->dot(*ykp1Rcp);
    double bTg = bRcp->dot(*ykp1Rcp);
    double gTAg = tempVecRcp->dot(*ykp1Rcp);
    std::cout << "             " << xTAg << " " << bTg << "  " << gTAg << std::endl;
    if (fabs(gTAg) < 1e-14) {
        gTAg += 1e-14; // prevent div 0 error
    }
    xkRcp->update(0.001 * (xTAg + 0.5 * bTg) / gTAg, *tempVecRcp, 1.0); //  xk=x0 - 0.001*(A.dot(xkm1)+b)
                                                                        // magic number 0.001 to start iteration

    ARcp->apply(*xkRcp, *ykRcp); // ykm1=A.dot(xkm1)
    mvCount++;
    MPI_Barrier(MPI_COMM_WORLD);
    ykmykm1Rcp->update(1.0, *ykRcp, -1.0, *ykm1Rcp, 0.0); // ykmykm1=yk-ykm1
    xkmxkm1Rcp->update(1.0, *xkRcp, -1.0, *xkm1Rcp, 0.0); // xkmxkm1=xk-xkm1

    int iteCount = 0;
    while (iteCount < iteMax) {
        iteCount++;
        MPI_Barrier(MPI_COMM_WORLD);

        // Barzilai-Borwein step size Choice 2
        // alphak=(xk-xkm1).dot(ykmykm1)/(ykmykm1).dot(ykmykm1)
        // double a = xkmxkm1Rcp->dot(*ykmykm1Rcp);
        // double b = ykmykm1Rcp->dot(*ykmykm1Rcp);

        // Barzilai-Borwein step size Choice 1, converges faster
        // alphak=(xk-xkm1).dot(xkmxkm1)/(xkmxkm1).dot(ykmykm1)
        double a = xkmxkm1Rcp->dot(*xkmxkm1Rcp);
        double b = xkmxkm1Rcp->dot(*ykmykm1Rcp);

        if (fabs(b) < 1e-14) {
            b = b + 1e-14; // prevent div 0 error
        }
        double alphak = a / b;
        // for debug
        // std::cout << a << " " << b << " " << alphak << std::endl;
        if (alphak < 1e-9) {
            // step size too small
            history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, alphak, 0}});
            break;
        }

        // xkp1=xk-alphak*(yk+b) # gradient descent
        xkp1Rcp->update(alphak, *ykRcp, alphak, *bRcp, 0.0);
        xkp1Rcp->update(1.0, *xkRcp, -1.0);
        clipZero(xkp1Rcp); // Projection xkp1 >= 0

        // update
        ARcp->apply(*xkp1Rcp, *ykp1Rcp); // ykp1=A.dot(xkp1)
        mvCount++;
        ykmykm1Rcp->update(1.0, *ykp1Rcp, -1.0, *ykRcp, 0.0); // yerror
        xkmxkm1Rcp->update(1.0, *xkp1Rcp, -1.0, *xkRcp, 0.0); // xerror
#ifdef DEBUGCOLLCP
        // detailed tolerance for debug build
        // check convergence
        double resxMax = xkmxkm1Rcp->normInf();
        double resAxbMax = ykmykm1Rcp->normInf();
        double resPhi = checkResiduePhi(xkp1Rcp, ykp1Rcp, bRcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, resxMax, resAxbMax, alphak, resPhi, 1.0 * mvCount}});
        if (fabs(resxMax) < tol && fabs(resAxbMax) < tol && fabs(AxbDotx) < tol) {
            break;
        }
#else
        // use simple phi tolerance check
        double resPhi = checkResiduePhi(xkp1Rcp, ykp1Rcp, bRcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, alphak, resPhi, 1.0 * mvCount}});
        if (fabs(resPhi) < tol) {
            break;
        }
#endif

        // prepare next iteration
        // swap the contents of pointers directly, be careful
        ykm1Rcp.swap(ykRcp); // ykm1=yk
        xkm1Rcp.swap(xkRcp); // xkm1=xk
        ykRcp.swap(ykp1Rcp); // yk=ykp1, ykp1 to be updated;
        xkRcp.swap(xkp1Rcp); // xk=xkp1, xkp1 to be updated;

        // ykm1Rcp->scale(1.0, *ykRcp);  // ykm1=yk
        // xkm1Rcp->scale(1.0, *xkRcp);  // xkm1=xk

        // ykRcp->scale(1.0, *ykp1Rcp);  // yk=ykp1
        // xkRcp->scale(1.0, *xkp1Rcp);  // xk=xkp1
    }
    xsolRcp = xkp1Rcp; // return solution

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

int CPSolver::LCP_APGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const {
    MPI_Barrier(MPI_COMM_WORLD);
    int mvCount = 0;
    if (myRank == 0) {
        std::cout << "solving APGD" << std::endl;
        std::cout << "ARcp" << ARcp->description() << std::endl;
    }
    // map must match
    TEUCHOS_TEST_FOR_EXCEPTION(!this->mapRcp->isSameAs(*(xsolRcp->getMap())), std::invalid_argument,
                               "xsolrcp and A operator do not have the same Map.");
    Teuchos::RCP<TV> xkRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy
    Teuchos::RCP<TV> ykRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy, yk=xk
    Teuchos::RCP<TV> xhatkRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    xhatkRcp->putScalar(1.0);

    double thetak = 1;
    double thetakp1 = 1;
    Teuchos::RCP<TV> xkdiffRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    xkdiffRcp->update(1.0, *(xhatkRcp.getConst()), -1.0, *(xkRcp.getConst()), 0.0);
    std::cout << "initial vector allocated" << std::endl;

    Teuchos::RCP<TV> tempVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));

    // std::cout << "tempVec allocated" << std::endl;
    ARcp->apply(*xkdiffRcp, *tempVecRcp);
    mvCount++;
    // std::cout << "A op applied" << std::endl;

    const double tempNorm2 = tempVecRcp->norm2();
    const double xkdiffNorm2 = xkdiffRcp->norm2();
    double Lk = (tempNorm2 / xkdiffNorm2);
    double tk = 1.0 / Lk;
    // std::cout << "initial step length calculated" << std::endl;
    // std::cout << "tk: " << tk << std::endl;

    // allocate other vectors
    Teuchos::RCP<TV> gVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> xkp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> ykp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> AxbRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> Axbkp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> resxRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> resAxbRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    history.push_back(std::array<double, 6>{{0, 0, 0, tk, 0, 1.0 * mvCount}});

    ARcp->apply(*xsolRcp, *AxbRcp); // Ax
    mvCount++;
    AxbRcp->update(1.0, *bRcp, 1.0); // Ax+b
    // enter main loop
    int iteCount = 0;
    while (iteCount < iteMax) {
        iteCount++;
        MPI_Barrier(MPI_COMM_WORLD);
        ARcp->apply(*ykRcp, *gVecRcp);
        mvCount++;
        gVecRcp->update(1.0, *bRcp, 1.0); // g=A.dot(yk)+b
        bool smallStepFlag = false;
        while (smallStepFlag == false) {
            // update and projection, xkp1=(yk-tk*g).clip(min=0)
            xkp1Rcp->update(1.0, *ykRcp, -tk, *gVecRcp, 0);
            clipZero(xkp1Rcp);
            // dumpTV(xkp1Rcp, "xkp1Proj");

            //  xkdiff=xkp1-yk
            xkdiffRcp->update(1.0, *xkp1Rcp, -1.0, *ykRcp, 0.0);

            //  calc Lifshitz condition
            ARcp->apply(*xkp1Rcp, *tempVecRcp);
            mvCount++;
            Axbkp1Rcp->update(1.0, *tempVecRcp, 1.0, *bRcp, 0.0); // Ax+b
            double leftTerm1 = xkp1Rcp->dot(*tempVecRcp) * 0.5;   // xkp1.dot(A.dot(xkp1))*0.5
            double leftTerm2 = xkp1Rcp->dot(*bRcp);               // xkp1.dot(b)

            ARcp->apply(*ykRcp, *tempVecRcp);
            mvCount++;
            double rightTerm1 = ykRcp->dot(*tempVecRcp) * 0.5;         // yk.dot(A.dot(yk))*0.5
            double rightTerm2 = ykRcp->dot(*bRcp);                     // yk.dot(b)
            double rightTerm3 = gVecRcp->dot(*xkdiffRcp);              // g.dot(xkdiff)
            double rightTerm4 = 0.5 * Lk * pow(xkdiffRcp->norm2(), 2); // 0.5*Lk*(xkdiff).dot(xkdiff)
            if ((leftTerm1 + leftTerm2) < (rightTerm1 + rightTerm2 + rightTerm3 + rightTerm4)) {
                break;
            }
            Lk *= 2;
            tk = 1 / Lk;
            if (tk < 1e-8) {
                smallStepFlag = true;
            }
        }
        if (smallStepFlag == true) {
            history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, tk, 0, 1.0 * mvCount}});
            break; // 1 for message 'step too small'
        }

        // APGD iterations
        resxRcp->update(1.0, *xkp1Rcp, -1.0, *xkRcp, 0.0);      // resx=xkp1-xk
        resAxbRcp->update(1.0, *Axbkp1Rcp, -1.0, *AxbRcp, 0.0); // resAxb = A.dot(xkp1) - A.dot(xkp)
        double AxbDotx = Axbkp1Rcp->dot(*xkp1Rcp);
#ifdef DEBUGCOLLCP
        // detailed tolerance for debug build
        // check convergence
        double resxMax = resxRcp->normInf();
        double resAxbMax = resAxbRcp->normInf();
        double resPhi = checkResiduePhi(xkp1Rcp, Axbkp1Rcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, resxMax, resAxbMax, tk, resPhi, 1.0 * mvCount}});
        if (fabs(resxMax) < tol && fabs(resAxbMax) < tol && fabs(resPhi) < tol) {
            break;
        }
#else
        // use simple phi tolerance check
        double resPhi = checkResiduePhi(xkp1Rcp, Axbkp1Rcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, tk, resPhi, 1.0 * mvCount}});
        if (fabs(resPhi) < tol) {
            break;
        }
#endif

        if (gVecRcp->dot(*resxRcp) > 0) {
            ykp1Rcp->scale(1.0, *xkp1Rcp); // ykp1=xkp1
            thetakp1 = 1;
        } else {
            thetakp1 = (-thetak * thetak + thetak * sqrt(4 + thetak * thetak)) / 2;
            double betakp1 = thetak * (1 - thetak) / (thetak * thetak + thetakp1);
            // ykp1=xkp1+betakp1*(xkp1-xk)
            ykp1Rcp->update((1 + betakp1), *xkp1Rcp, -betakp1, *xkRcp, 0);
        }

        // next iteration
        // swap the contents of pointers directly, be careful
        ykRcp.swap(ykp1Rcp);    // yk=ykp1, ykp1 to be updated;
        xkRcp.swap(xkp1Rcp);    // xk=xkp1, xkp1 to be updated;
        AxbRcp.swap(Axbkp1Rcp); // Axbk=Axbkp1, Axbkp1 to be updated;

        thetak = thetakp1;
        Lk *= 0.9;
        tk = 1 / Lk;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    xsolRcp = xkp1Rcp;

    return 0;
}

// class diagOperator : public TOP {
//  private:
//   Teuchos::RCP<TV> diagRcp;

//  public:
//   diagOperator(const Teuchos::RCP<const TOP> &AopRcp) {
//     // allocate, fill with identity
//     diagRcp = Teuchos::rcp(new TV(AopRcp->getDomainMap(), false));
//     diagRcp.putScalar(1);
//   }
// }

class mmNewtonOperator : public TOP {
  private:
    Teuchos::RCP<const TOP> AopRcp;
    Teuchos::RCP<const TV> maskRcp;
    Teuchos::RCP<TCMAT> precOpRcp;

  public:
    mmNewtonOperator(const Teuchos::RCP<const TOP> &AopRcp_, const Teuchos::RCP<const TV> &maskRcp_) : AopRcp(AopRcp_) {
        maskUpdate(maskRcp_);
    }

    void maskUpdate(const Teuchos::RCP<const TV> &maskRcp_) {
        maskRcp = maskRcp_;
        // update precOp
    }

    void opUpdate(const Teuchos::RCP<const TOP> &AopRcp_) { AopRcp = AopRcp_; }

    ~mmNewtonOperator(){};

    Teuchos::RCP<const TMAP> getDomainMap() const {
        return this->AopRcp->getDomainMap(); // Get the domain Map of this Operator subclass.
    }
    Teuchos::RCP<const TMAP> getRangeMap() const {
        return this->AopRcp->getRangeMap(); // Get the range Map of this Operator subclass.
    }

    bool hasTransposeApply() const { return false; }

    // Compute Y := alpha Op X + beta Y.
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const {
        // step 1 Y=A*x
        AopRcp->apply(X, Y);
        const size_t numVecsX = X.getNumVectors();
        const size_t numVecsY = Y.getNumVectors();
        assert(numVecsX == numVecsY);
        const int localSize = X.getMap()->getNodeNumElements();
        // step 2  apply the mask

        auto xView = X.getLocalView<Kokkos::HostSpace>();
        auto yView = Y.getLocalView<Kokkos::HostSpace>();
        auto maskView = maskRcp->getLocalView<Kokkos::HostSpace>();
        Y.modify<Kokkos::HostSpace>();
        for (int c = 0; c < xView.dimension_1(); c++) {
            for (int i = 0; i < xView.dimension_0(); i++) {
                if (maskView(i, c) > 0.5) {
                    yView(i, c) = xView(i, c);
                }
            }
        }
    }
};

int CPSolver::LCP_mmNewton(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const {
    MPI_Barrier(MPI_COMM_WORLD);
    int mvCount = 0;
    if (myRank == 0) {
        std::cout << "solving mmNewton" << std::endl;
        std::cout << "ARcp" << ARcp->description() << std::endl;
    }
    // map must match
    TEUCHOS_TEST_FOR_EXCEPTION(!this->mapRcp->isSameAs(*(xsolRcp->getMap())), std::invalid_argument,
                               "xsolrcp and A operator do not have the same Map.");

    Teuchos::RCP<TV> xRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy, xk=x0
    Teuchos::RCP<TV> yRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> xkRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> ykRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));

    Teuchos::RCP<TV> tempVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));

    Teuchos::RCP<TV> dxRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));      // zero out
    Teuchos::RCP<TV> nablaHRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false)); // zero out

    Teuchos::RCP<TV> HmmRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));   // minimum map
    Teuchos::RCP<TV> HmaskRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false)); // minimum map

    // Magic constants
    const double alpha = 0.5;
    const double beta = 0.001;
    const double gamma = 1e-28;
    const double eps = 1e-14;
    const double rho = 1e-14;
    const double gmres_tol = tol;
    const double tol_abs = 1e-14;
    const double tol_rel = tol;
    double err = 1e20;
    int iteCount = 0;

    // allocate GMRES space
    Teuchos::RCP<mmNewtonOperator> maskAOpRcp = Teuchos::rcp(new mmNewtonOperator(ARcp, HmaskRcp));
    // set Belos object
    Belos::SolverFactory<TOP::scalar_type, TMV, TOP> factory;
    // Make an empty new parameter list.
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Num Blocks", 50); // usless for BicgStab. the restart m for GMRES.
    solverParams->set("Maximum Iterations", 100);
    solverParams->set("Convergence Tolerance", tol * 10);
#ifdef DEBUGCOLLCP
    solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);
#else
    solverParams->set("Verbosity", Belos::Errors + Belos::Warnings);
#endif
    auto solverRCP = factory.create("GMRES", solverParams); // Create the GMRES solver.
    auto problemRCP = Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(maskAOpRcp, dxRcp, HmmRcp));
    problemRCP->setProblem(); // necessary
    solverRCP->setProblem(problemRCP);

    // TODO: Preconditioner
    // Teuchos::RCP<TOP> precOp;
    // problemRCP->setRightPrec(precOp);
    ARcp->apply(*xRcp, *yRcp);
    mvCount++;
    yRcp->update(1.0, *bRcp, 1.0); // y = A.dot(x) + b
    double oldErr;
    double tk = 0;
    while (iteCount++ < iteMax) {
    // check convergence
    // xkRcp, ykRcp saves last step

#ifdef DEBUGCOLLCP
        // detailed tolerance for debug build
        // check convergence
        tempVecRcp->update(1.0, *xRcp, -1.0, *xkRcp, 0.0);
        double resxMax = tempVecRcp->normInf();
        tempVecRcp->update(1.0, *yRcp, -1.0, *ykRcp, 0.0);
        double resAxbMax = tempVecRcp->normInf();
        double resPhi = checkResiduePhi(xkRcp, ykRcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, resxMax, resAxbMax, tk, resPhi, 1.0 * mvCount}});
        if (tk > 0 && fabs(resxMax) < tol && fabs(resAxbMax) < tol && fabs(resPhi) < tol) {
            // converge, stop
            break;
        }
#else
        // use simple phi tolerance check
        double resPhi = checkResiduePhi(xkRcp, ykRcp, tempVecRcp);
        history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, tk, resPhi, 1.0 * mvCount}});
        if (tk > 0 && fabs(resPhi) < tol) {
            break;
        }
#endif

        // minimum map, H = minmap(x,y)
        hMinMap(xRcp, yRcp, HmmRcp, HmaskRcp);
        oldErr = err; // old_err = err
        // Calculate merit value, error: err = 0.5*H.dot(H)
        err = 0.5 * pow(HmmRcp->norm2(), 2);

        //##### Test the stopping criterias used #####
        double rel_err = fabs(err - oldErr) / fabs(oldErr);
        if (rel_err < tol_rel || err < tol_abs) {
            break;
        }

        /*
                ##### Solving the Newton system
            restart = min(N, 20) # Number of iterates done before Restart
                                 # for GMRES should restart
            S = np.where(y < x)
            J = np.identity(N)
            J[S,:] = A[S,:]
            dx = np.zeros(N)
            dx = gmres(J, (-H), tol=gmres_tol, restart=restart)[0].reshape(N)
        */
        maskAOpRcp->maskUpdate(HmaskRcp.getConst());
        // dxRcp->putScalar(0);
        problemRCP->setProblem(); // necessary to update the solver
        // solverRCP->setProblem(problemRCP);
        MPI_Barrier(MPI_COMM_WORLD);
        Belos::ReturnType result = solverRCP->solve(); //  maskAOpRcp *dx = -H, solution is at dxRcp
        mvCount += solverRCP->getNumIters();
        MPI_Barrier(MPI_COMM_WORLD);
        dxRcp->scale(-1.0);

        // nabla_H = H.dot(J) # descent direction
        maskAOpRcp->apply(*HmmRcp, *nablaHRcp);

        //    # Tests whether the search direction is below machine precision.
        if (dxRcp->normInf() < tol_abs) {
            std::cout << "Search direction below machine precision at iteration " << iteCount << std::endl;
            std::cout << "Using descent direction" << std::endl;
            dxRcp->update(-1.0, *nablaHRcp, 0.0); // dx = -nabla_H
        }

        // # Test whether we are stuck in a local minima
        if (nablaHRcp->norm2() < tol_abs) {
            std::cout << "local minimum" << std::endl;
            return 1;
        }

        //    # Test whether our direction is a sufficient descent direction
        if (nablaHRcp->dot(*dxRcp) > -rho * pow(dxRcp->norm2(), 2)) {
            std::cout << "Non descend direction at iteration " << iteCount << ", choosing gradient as search direction"
                      << std::endl;
            dxRcp->update(-1.0, *nablaHRcp, 0.0);
        }

        // ##### Armijo backtracking combined with a projected line-search #####
        double tau = 1.0;
        double f0 = err;
        double gradf = beta * nablaHRcp->dot(*dxRcp);
        // # Perform backtracking line search
        while (1) {
            xkRcp->update(1.0, *xRcp, tau, *dxRcp, 0.0); // xk=x + dx*tau
            clipZero(xkRcp);                             // xk=max(0,xkRcp);
            ARcp->apply(*xkRcp, *ykRcp);
            mvCount++;
            ykRcp->update(1.0, *bRcp, 1.0);            // yk = y_k = np.dot(A,x_k)+b
            hMinMap(ykRcp, xkRcp, HmmRcp, HmaskRcp);   // H_k = minmap(y_k,x_k)
            double fk = 0.5 * pow(HmmRcp->norm2(), 2); // f_k = 0.5*(H_k.dot(H_k))

            if (fk <= f0 + tau * gradf) {
                // # Test Armijo condition for sufficient decrease
                break;
            }
            // # Test whether the stepsize has become too small
            if (tau * tau < gamma) {
                break;
            }
            tau *= alpha;
        }
        // # Update iterate with result from line search.
        tk = tau;
        xRcp.swap(xkRcp);
        yRcp.swap(ykRcp);
    }

    xsolRcp.swap(xRcp);

    return 0;
}

int CPSolver::test_LCP(double tol, int maxIte, int solverChoice) {
    IteHistory history;

    Teuchos::RCP<TV> xsolRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // zero initial guess

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    switch (solverChoice) {
    case 0:
        LCP_mmNewton(xsolRcp, tol, maxIte, history);
        break;
    case 1:
        LCP_APGD(xsolRcp, tol, maxIte, history);
        break;
    case 2:
        LCP_BBPGD(xsolRcp, tol, maxIte, history);
        break;
    case 3:
        LCP_APGD(xsolRcp, 100 * tol, maxIte, history);
        LCP_mmNewton(xsolRcp, tol, maxIte, history);
        break;
    case 4:
        LCP_BBPGD(xsolRcp, 100 * tol, maxIte, history);
        LCP_mmNewton(xsolRcp, tol, maxIte, history);
        break;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Solving time: " << time_span.count() << " seconds." << std::endl;

    // test result
    Teuchos::RCP<TV> AxbRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    ARcp->apply(*xsolRcp, *AxbRcp);  // Ax
    AxbRcp->update(1.0, *bRcp, 1.0); // Ax+b
    dumpTV(xsolRcp, "xsol.mtx");
    dumpTV(AxbRcp, "Axb.mtx");
    auto xView = xsolRcp->getLocalView<Kokkos::HostSpace>();
    auto yView = AxbRcp->getLocalView<Kokkos::HostSpace>();

    int xerrorN = 0;
    double xerrormax = 0.0;
    int yerrorN = 0;
    double yerrormax = 0.0;
    assert(xView.dimension_0() == localSize);
    for (int c = 0; c < xView.dimension_1(); c++) {
        for (int i = 0; i < xView.dimension_0(); i++) {
            if (xView(i, c) < 0) {
                xerrorN++;
                xerrormax = std::max(xerrormax, static_cast<double>(-xView(i, c)));
            }
            if (yView(i, c) < 0) {
                yerrorN++;
                yerrormax = std::max(yerrormax, static_cast<double>(-yView(i, c)));
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myRank == 0) {
        // std::cout << myRank << std::endl;
        MPI_Reduce(MPI_IN_PLACE, &xerrorN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &yerrorN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &xerrormax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &yerrormax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // std::cout << myRank << std::endl;
    } else {
        // std::cout << myRank << std::endl;
        int tempint;
        double tempdouble;
        MPI_Reduce(&xerrorN, &tempint, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&yerrorN, &tempint, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&xerrormax, &tempdouble, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&yerrormax, &tempdouble, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        // std::cout << myRank << std::endl;
    }

    if (myRank == 0) {
        for (auto &p : history) {
            std::cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << p[5] << std::endl;
        }
    }

    if (myRank == 0) {
        std::cout << "solution quality: " << std::endl;
        std::cout << "negative x N: " << std::scientific << xerrorN << std::endl;
        std::cout << "xerror max: " << std::scientific << xerrormax << std::endl;
        std::cout << "negative y N: " << std::scientific << yerrorN << std::endl;
        std::cout << "yerror max: " << std::scientific << yerrormax << std::endl;
        std::cout << "Solving time: " << time_span.count() << " seconds." << std::endl;
        // double dotyTx = AxbRcp->dot(*xsolRcp);
        // std::cout << "y^T x: " << std::scientific << dotyTx << std::endl;
        // stuck for >=2 processors, unknown
    }

    return 0;
}
