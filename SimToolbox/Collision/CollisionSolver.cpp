#include "CollisionSolver.hpp"

// The number of collision blocks can be zero on one, several or all nodes
// be careful about this special case.

void CollisionSolver::setup(CollisionBlockPool &collision_, Teuchos::RCP<TMAP> &objMobMapRcp_, double dt_,
                            double bufferGap_) {
    reset();
    objMobMapRcp = objMobMapRcp_;
    forceColRcp = Teuchos::rcp(new TV(objMobMapRcp, true));
    velocityColRcp = Teuchos::rcp(new TV(objMobMapRcp, true));

    // check global empty
    // bool emptyFlagLocal = collisionLocalIsEmpty(collision_);
    // bool emptyFlagGlobal = true;
    // MPI_Allreduce(&emptyFlagLocal, &emptyFlagGlobal, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
    // if (emptyFlagGlobal) {
    //     // TMAP seems does not allow 0 global size. Add a manual fix
    // }

    setupCollisionBlockQueThreadIndex(collision_);

    // step 1 setup maps
    auto commRcp = objMobMapRcp->getComm();
    gammaMapRcp = getTMAPFromLocalSize(queueThreadIndex.back(), commRcp);

    // const int nObjLocal = objMobMapRcp->getNodeNumElements() / 6;
    // const int nObjGlobal = objMobMapRcp->getGlobalNumElements() / 6;
    // TEUCHOS_ASSERT(nObjLocal * 6 == objMobMapRcp->getNodeNumElements());
    // TEUCHOS_ASSERT(nObjGlobal * 6 == objMobMapRcp->getGlobalNumElements());

    if (gammaMapRcp->getGlobalNumElements() == 0) {
        if (commRcp->getRank() == 0)
            std::cout << "No collision detected, will set forceCol and velCol to zero" << std::endl;
    } else {
        // setup vecs
        gammaRcp = Teuchos::rcp(new TV(gammaMapRcp, true));

        // step 2 setup FcTrans
        // assert(velocity.size() == nObjLocal * 6);
        // setupVnVec(collision, velocity);
        setupFcTrans(collision_);

        // step 3 setup the const b in LCP, stored in phi0. need FcTrans and Vn
        setupPhi0Vec(collision_, dt_, bufferGap_);
    }
}

void CollisionSolver::solveCollision(Teuchos::RCP<TOP> &matMobilityRcp_, Teuchos::RCP<TV> &velocityKnownRcp_) {
    // should return valid (0) results for no collisions
    this->matMobilityRcp = matMobilityRcp_;
    this->vnRcp = velocityKnownRcp_;

    if (gammaMapRcp->getGlobalNumElements() == 0) {
        // global no collision.
        forceColRcp->putScalar(0);
        velocityColRcp->putScalar(0);
        return;
    }

    // set the const b
    setupBVec();

    // create the solver
    IteHistory history;
    // the operator A in LCP
    Teuchos::RCP<CPMatOp> AmatRcp = Teuchos::rcp(new CPMatOp(matMobilityRcp, matFcTransRcp));
    auto commRcp = objMobMapRcp->getComm();
    CPSolver myLCPSolver(AmatRcp, bRcp);
    if (commRcp->getRank() == 0) {
        printf("start solving\n");
    }
    myLCPSolver.LCP_BBPGD(gammaRcp, res, maxIte, history);

    if (commRcp->getRank() == 0 && history.size() > 0) {
        auto &p = history.back();
        std::cout << "LCP residue: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << std::endl;
    }

    if (newton) {
        myLCPSolver.LCP_mmNewton(gammaRcp, res * 0.01, maxIte, history);
    }

    // save the solution
    forceColRcp = AmatRcp->forceVecRcp;
    velocityColRcp = AmatRcp->velVecRcp;
    if (commRcp->getRank() == 0)
        printf("col force velocity saved\n");

#ifdef DEBUGLCPCOL
    if (commRcp->getRank() == 0 && history.size() > 0) {
        for (auto &p : history) {
            std::cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << std::endl;
        }
    }
    dumpTV(forceColRcp, "forceCol");
    dumpTV(gammaRcp, "gammaSolLCP");
#endif
}

// general FcTrans mat: 6 mobility dof per object to 1 dof per collision
// each row: 6 nnz if collision with boundary, 12 nnz if collision between objects
// each 6nnz for an object: gx.ux+gy.uy+gz.uz+(gzpy-gypz)wx+(gxpz-gzpx)wy+(gypx-gxpy)wz
// FOR SIMPLICITY of counting nnz, always use 12 nnz per collision. Fill 0s when only 6 nnz is necessary
void CollisionSolver::setupFcTrans(CollisionBlockPool &collision_) {

#ifdef DEBUGLCPCOL
    dumpCollision(collision_);
#endif

    const int nThreads = collision_.size(); // each thread process one queue in the pool
    const int localGammaSize = queueThreadIndex.back();

    // initialize the Fc^T collision matrix
    // rowMap = colBlockMap, colMap=mobMatMap
    // six entries for each row for spheres (each colBlock), gI[x,y,z], gJ[x,y,z]
    Kokkos::View<size_t *> rowPointers("rowPointers", localGammaSize + 1);
    rowPointers[0] = 0;
    // collision blocks
    int rowPointerIndex = 0;
    std::vector<int> colIndexThread(collision_.size());
    for (int i = 0; i < collision_.size(); i++) {
        const int jsize = collision_[i].size();
        int colIndexCount = 0;
        for (int j = 0; j < jsize; j++) {
            rowPointerIndex++;
            const int colBlockNNZ = (collision_[i][j].oneSide ? 6 : 12);
            rowPointers[rowPointerIndex] = rowPointers[rowPointerIndex - 1] + colBlockNNZ;
            colIndexCount += colBlockNNZ;
        }
        colIndexThread[i] = colIndexCount;
    }
    if (rowPointerIndex != localGammaSize) {
        printf("rowPointerIndexError in collision solver\n");
        exit(1);
    }

    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localGammaSize]);
    Kokkos::View<double *> values("values", rowPointers[localGammaSize]);

#pragma omp parallel for num_threads(nThreads)
    for (int threadId = 0; threadId < nThreads; threadId++) {
        // each thread process a queue
        const auto &colBlockQue = collision_[threadId];
        const int colBlockNum = colBlockQue.size();
        const int colBlockIndexBase = queueThreadIndex[threadId];
        // obj-obj collision: 12 nnz
        // obj-boundary collision: 6 nnz
        // find proper kk
        int kk = 0;
        for (int t = 0; t < threadId; t++) {
            kk += colIndexThread[t];
        }

        for (int j = 0; j < colBlockNum; j++) {
            columnIndices[kk] = 6 * colBlockQue[j].globalIndexI;
            columnIndices[kk + 1] = 6 * colBlockQue[j].globalIndexI + 1;
            columnIndices[kk + 2] = 6 * colBlockQue[j].globalIndexI + 2;
            columnIndices[kk + 3] = 6 * colBlockQue[j].globalIndexI + 3;
            columnIndices[kk + 4] = 6 * colBlockQue[j].globalIndexI + 4;
            columnIndices[kk + 5] = 6 * colBlockQue[j].globalIndexI + 5;
            // each 6nnz for an object: gx.ux+gy.uy+gz.uz+(gzpy-gypz)wx+(gxpz-gzpx)wy+(gypx-gxpy)wz
            // 6 nnz for I
            {
                const double &gx = colBlockQue[j].normI[0];
                const double &gy = colBlockQue[j].normI[1];
                const double &gz = colBlockQue[j].normI[2];
                const double &px = colBlockQue[j].posI[0];
                const double &py = colBlockQue[j].posI[1];
                const double &pz = colBlockQue[j].posI[2];
                values[kk] = gx;
                values[kk + 1] = gy;
                values[kk + 2] = gz;
                values[kk + 3] = (gz * py - gy * pz);
                values[kk + 4] = (gx * pz - gz * px);
                values[kk + 5] = (gy * px - gx * py);
            }
            if (!colBlockQue[j].oneSide) {
                columnIndices[kk + 6] = 6 * colBlockQue[j].globalIndexJ;
                columnIndices[kk + 7] = 6 * colBlockQue[j].globalIndexJ + 1;
                columnIndices[kk + 8] = 6 * colBlockQue[j].globalIndexJ + 2;
                columnIndices[kk + 9] = 6 * colBlockQue[j].globalIndexJ + 3;
                columnIndices[kk + 10] = 6 * colBlockQue[j].globalIndexJ + 4;
                columnIndices[kk + 11] = 6 * colBlockQue[j].globalIndexJ + 5;
                // 6 nnz for J, should be 0 if is a boundary. set normJ and posJ to zero before calling this
                {
                    const double &gx = colBlockQue[j].normJ[0];
                    const double &gy = colBlockQue[j].normJ[1];
                    const double &gz = colBlockQue[j].normJ[2];
                    const double &px = colBlockQue[j].posJ[0];
                    const double &py = colBlockQue[j].posJ[1];
                    const double &pz = colBlockQue[j].posJ[2];
                    values[kk + 6] = gx;
                    values[kk + 7] = gy;
                    values[kk + 8] = gz;
                    values[kk + 9] = (gz * py - gy * pz);
                    values[kk + 10] = (gx * pz - gz * px);
                    values[kk + 11] = (gy * px - gx * py);
                }
                kk += 12;
            } else {
                kk += 6;
            }
        }
    }

    auto commRcp = objMobMapRcp->getComm();
    // Teuchos::RCP<TMAP> colMapRcp = getFullCopyTMAPFromGlobalSize(objMobMapRcp->getGlobalNumElements(), commRcp);
    const int colIndexCount = rowPointers[localGammaSize];
    std::vector<int> colMapIndex(colIndexCount);
#pragma omp parallel for
    for (int i = 0; i < colIndexCount; i++) {
        colMapIndex[i] = columnIndices[i];
    }

    // sort and unique
    std::sort(colMapIndex.begin(), colMapIndex.end());
    colMapIndex.erase(std::unique(colMapIndex.begin(), colMapIndex.end()), colMapIndex.end());

    Teuchos::RCP<TMAP> colMapRcp = Teuchos::rcp(
        new TMAP(Teuchos::OrdinalTraits<int>::invalid(), colMapIndex.data(), colMapIndex.size(), 0, commRcp));

    // convert columnIndices from global column index to local column index according to colMap
    auto &colmap = *colMapRcp;
#pragma omp parallel for
    for (int i = 0; i < colIndexCount; i++) {
        columnIndices[i] = colmap.getLocalElement(columnIndices[i]);
    }

    matFcTransRcp = Teuchos::rcp(new TCMAT(gammaMapRcp, colMapRcp, rowPointers, columnIndices, values));
    matFcTransRcp->fillComplete(objMobMapRcp, gammaMapRcp); // domainMap, rangeMap

#ifdef DEBUGLCPCOL
    std::cout << "FcTransConstructed: " << matFcTransRcp->description() << std::endl;
    dumpTCMAT(matFcTransRcp, "matFcTrans");
#endif

    return;
}

void CollisionSolver::setupCollisionBlockQueThreadIndex(CollisionBlockPool &collision_) {
    const int poolSize = collision_.size();
    queueThreadIndex.resize(poolSize + 1, 0);
    for (int i = 1; i <= poolSize; i++) {
        queueThreadIndex[i] = queueThreadIndex[i - 1] + collision_[i - 1].size();
    }
}

void CollisionSolver::dumpCollision(CollisionBlockPool &collision_) const {
#ifdef DEBUGLCPCOL
    // dump collision data
    for (auto &colBlockQue : collision_) {
        std::cout << colBlockQue.size() << " collisions in this queue" << std::endl;
        for (auto &colBlock : colBlockQue) {
            std::cout << colBlock.globalIndexI << " " << colBlock.globalIndexJ << "  sep:" << colBlock.phi0
                      << std::endl;
        }
    }
#endif
}

// current known value of constraints
void CollisionSolver::setupPhi0Vec(CollisionBlockPool &collision_, double dt_, double bufferGap_) {

    phi0Rcp = Teuchos::rcp(new TV(gammaMapRcp.getConst(), false));

    // raw ptrs to local values
    auto phi0_2d = phi0Rcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    phi0Rcp->modify<Kokkos::HostSpace>();

    const int nThreads = collision_.size();

#pragma omp parallel for num_threads(nThreads)
    for (int threadId = 0; threadId < nThreads; threadId++) {
        // each thread process a queue
        const auto &colBlockQue = collision_[threadId];
        const int colBlockNum = colBlockQue.size();
        const int colBlockIndexBase = queueThreadIndex[threadId];
        for (int j = 0; j < colBlockNum; j++) {
            const double phi0buf = colBlockQue[j].phi0 / dt_ + bufferGap_;
            phi0_2d(colBlockIndexBase + j, 0) = phi0buf;
        }
    }

#ifdef DEBUGLCPCOL
    std::cout << "phi0 vector Constructed: " << phi0Rcp->description() << std::endl;
    dumpTV(phi0Rcp, "phi0Vec");
#endif
}

// initial guess of unknown gamma
void CollisionSolver::setupGammaVec(CollisionBlockPool &collision_) {
    gammaRcp = Teuchos::rcp(new TV(gammaMapRcp.getConst(), false));
    // raw ptrs to local values
    auto gamma_2d = gammaRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    gammaRcp->modify<Kokkos::HostSpace>();

    const int nThreads = collision_.size();

#pragma omp parallel for num_threads(nThreads)
    for (int threadId = 0; threadId < nThreads; threadId++) {
        // each thread process a queue
        const auto &colBlockQue = collision_[threadId];
        const int colBlockNum = colBlockQue.size();
        const int colBlockIndexBase = queueThreadIndex[threadId];
        for (int j = 0; j < colBlockNum; j++) {
            gamma_2d(colBlockIndexBase + j, 0) = colBlockQue[j].gamma;
        }
    }

#ifdef DEBUGLCPCOL
    std::cout << "gamma vector Constructed: " << gammaRcp->description() << std::endl;
    dumpTV(gammaRcp, "gammaVec");
#endif
}

// known velocity
void CollisionSolver::setupVnVec(CollisionBlockPool &collision_, std::vector<double> &velocity_) {
    auto commRcp = objMobMapRcp->getComm();
    vnRcp = getTVFromVector(velocity_, commRcp);
#ifdef DEBUGLCPCOL
    std::cout << "Vn vector Constructed: " << vnRcp->description() << std::endl;
    dumpTV(vnRcp, "VnVec");
#endif
}

int CollisionSolver::getNumNegElements(Teuchos::RCP<TV> &vecRcp_) const {
    // count the total number of negtive phi0.
    int negNumber = 0;
    auto vec_2d = vecRcp_->getLocalView<Kokkos::HostSpace>(); // LeftLayout
#pragma omp parallel for reduction(+ : negNumber)
    for (int j = 0; j < vec_2d.dimension_0(); j++) {
        negNumber += (vec_2d(j, 0) < 0 ? 1 : 0);
    }
    // MPI_Allreduce(MPI_IN_PLACE, &negNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // return negNumber;
    int negNumberAll = 0;
    Teuchos::reduceAll(*(objMobMapRcp->getComm()), Teuchos::REDUCE_SUM, 1, &negNumber, &negNumberAll);

    return negNumberAll;
}

void CollisionSolver::setupBVec() {
    bRcp = Teuchos::rcp(new TV(*phi0Rcp, Teuchos::Copy));
    // the vec b in LCP stored in phi0
    // b = phi0
    // b = b + Fc^T * Vn
    bRcp->getMap()->getComm()->barrier();
    matFcTransRcp->apply(*vnRcp, *bRcp, Teuchos::NO_TRANS, 1.0, 1.0);
}
