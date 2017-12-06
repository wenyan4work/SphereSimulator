
// std
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <string>
#include <utility>
// #include <iostream>

// namespace HydroSphere
#include "HydroSphere.h"
#include "HydroSphereType.hpp"

// helper routines
#include "Preconditioner.hpp"
#include "ZDD.hpp"

// for debug. NDEBUG defined in pvfmm.hpp
#undef NDEBUG
#include <cassert>

namespace HydroSphere {

template <class T>
inline MPI_Datatype GetDataType() {
    static MPI_Datatype type = MPI_DATATYPE_NULL;
    if (type == MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(sizeof(T), MPI_BYTE, &type);
        MPI_Type_commit(&type);
        // deprecated c++ binding
        // type = MPI::BYTE.Create_contiguous(sizeof(T));
        // type.Commit();
    }
    return type;
};

class HydroSphereOperator : public TOP {

  private:
    // This is an implementation detail; users don't need to see it.
    HydroSphereSystem *systemPtr;

  public:
    // Constructor
    explicit HydroSphereOperator(HydroSphereSystem *systemPtr_) : systemPtr{systemPtr_} {}

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
        : systemPtr{systemPtr_}, Tscale{Tscale_}, Lscale{Lscale_}, Escale{Escale_}, Fscale{Escale_ / Lscale_},
          Torscale{Escale_} {
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

// max N, periodic type, fmm parameter, iterative solver parameter
HydroSphereSystem::HydroSphereSystem(const int maxChebN, const FMM_Wrapper::PAXIS pbc, const int multi_order,
                                     const double tol, const Evec3 &boxLow_, const Evec3 &boxHigh_,
                                     const double stkReg_, const double KReg_)
    : myFMM(multi_order, 1000, 0, pbc), iteres(tol), stkReg(stkReg_), KReg(KReg_), totalChebN(0),
      commRcp(Tpetra::DefaultPlatform::getDefaultPlatform().getComm()), myRank(commRcp->getRank()),
      nProcs(commRcp->getSize()), sphereGidFindDD(1000) {
    // initialization order determined by order in declaration
    // Tpetra objects
    if (commRcp->getRank() != myRank || commRcp->getSize() != nProcs) {
        std::cout << " hydro mpi rank and size error " << std::endl;
        exit(1);
    }

    //	Kokkos::initialize();
    boxLow = boxLow_;
    boxHigh = boxHigh_;

    // initialize FMM object
    myFMM.FMM_SetBox(boxLow[0], boxHigh[0], boxLow[1], boxHigh[1], boxLow[2], boxHigh[2]);

    // Reserve space
    fdistValid = false; // whether the force density is valid.
    commRcp->barrier();
    std::cout << "HydroSystem initialized\n";
}

// HydroSphereSystem::~HydroSphereSystem() {
//
//  }

void HydroSphereSystem::addSphere(const RigidSphereIO &rigidSphereIO) {
    sphereIO.emplace_back(rigidSphereIO);
    sphere.emplace_back(rigidSphereIO);
    // set up sphere data structure
    auto &rigidSphere = sphere.back();

    rigidSphere.myIntegrator = getChebIntegrator(rigidSphere.sphereDataIO.meshDelta, rigidSphere.sphereDataIO.alpha);
    rigidSphere.nbBlocks.clear();
    rigidSphere.fdistIndexLocal = INVALIDINDEX; // index of the starting point of fdist in assembled linear system
    const int chebN = rigidSphere.myIntegrator->chebN;
    rigidSphere.fdist.resize(3 * chebN + 3);
    // rigidSphere.fdistlast.resize(3 * chebN + 3);
    // rigidSphere.fdistout.resize(3 * chebN + 3);
    // rigidSphere.k0fdist.resize(3 * chebN + 3);
    // rigidSphere.Vmdist.resize(3 * chebN + 3);
    // rigidSphere.K0Vmdist.resize(3 * chebN + 3);
    // rigidSphere.vdistself.resize(3 * chebN + 3);
    // rigidSphere.vdist.resize(3 * chebN + 3);
    // rigidSphere.vdistReg.resize(3 * chebN + 3);
    rigidSphere.clearAll();

    rigidSphere.selfBlock.resize(3 * chebN + 3, 3 * chebN + 3); // self Operator

    // must be called when the integrator pointer is valid.
    rigidSphere.calculateSelfBlock();
    rigidSphere.calculateMobilityMatrix();

    // postcondition: sphere self data is valid, sphere neighbor blocks are not valid

    return;
}

void HydroSphereSystem::syncSphereIn() {
    // copy sphereIO data to sphere, before solution
    // pre: valid data, exchange, unpack and sort have been called
    // post: sphere.fdataIO=sphereIO,
    const int sphereNumber = sphere.size();
    assert(sphereNumber == sphereIO.size());
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        assert(sphere[i].sphereDataIO.gid == sphereIO[i].gid);
        // assert(sphereIO[i].gid != INVALIDINDEX);
        sphereIO[i].xdot.setZero();
        sphereIO[i].tdot.setZero();
        sphere[i].sphereDataIO = sphereIO[i]; // copy force and
    }
    return;
}

void HydroSphereSystem::syncSphereInForceTorque() {
    // copy sphereIO force and torque to sphere
    const int sphereNumber = sphere.size();
    assert(sphereNumber == sphereIO.size());
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        assert(sphere[i].sphereDataIO.gid == sphereIO[i].gid);
        assert(sphereIO[i].gid != INVALIDINDEX);
        sphere[i].sphereDataIO.force = sphereIO[i].force;   // copy force and
        sphere[i].sphereDataIO.torque = sphereIO[i].torque; // copy force and
    }
    return;
}

void HydroSphereSystem::syncSphereOut() {
    // copy sphere data to sphereIO, after solution
    // pre: sphere xdot,tdot is valid
    // post: sphere xdot,tdot copied to sphereIO
    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        assert(sphere[i].sphereDataIO.gid == sphereIO[i].gid);
        assert(sphere[i].sphereDataIO.force == sphereIO[i].force);
        sphereIO[i].xdot = sphere[i].sphereDataIO.xdot;
        sphereIO[i].tdot = sphere[i].sphereDataIO.tdot;
    }
    return;
}

void HydroSphereSystem::getReadyForSolve(SolverType solverType) {
    commRcp->barrier();
    // normalize direction vector
    int sphereNumber = sphereIO.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        sphereIO[i].direction.normalize();
    }

    exchangeSphere();
    syncSphereIn();
    constructDiscretization();
    updatefdistIndex();
    updateFMMTree();
    if (solverType == SolverType::IMPLICIT) {
        updateNeighborBlocks(true); // update nbBlockMat and nbRegMat
        assembleFlowRegMatrix(0);
        assembleApproxImpSolverMatrix(this->iteres / 100);
    } else {
        updateNeighborBlocks(false);
        assembleFlowRegMatrix(this->iteres / 100);
    }

    // clear temporary space
    sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].clearTemporaryVec();
    }

    return;
}

void HydroSphereSystem::updateNeighborBlocks(bool updateNbBlockMat) {
    // pre: valid sphereDataIO.nbInfos, valid fdistIndex, valid fdistMapRcp
    // post: valid sphere.nbBlocks
    const int sphereNumber = sphere.size();
    std::vector<int> nbBufferIndex(sphereNumber + 1); // last = total nb number
    nbBufferIndex[0] = 0;
    for (int i = 1; i < sphereNumber + 1; i++) {
        nbBufferIndex[i] = nbBufferIndex[i - 1] + sphere[i - 1].sphereDataIO.nbInfos.size();
        // sphere.nbBlocks is not initialized here yet
    }

    const int nbNumber = nbBufferIndex.back();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        assert(sphere[i].fdistIndexLocal != INVALIDINDEX);
        sphere[i].nbBlocks.clear();
        for (int j = 0; j < sphere[i].sphereDataIO.nbInfos.size(); j++) {
            assert(sphere[i].sphereDataIO.nbInfos[j].gid != INVALIDINDEX);
            sphere[i].nbBlocks.emplace_back(sphere[i].sphereDataIO.nbInfos[j]);
        }
    }

    // step 1 update neighbor information according to nbinfo
    struct NbInfoFind {
        int chebN;
        int gid;
        int fdistIndexGlobal;
        double alpha;
        double direction[3];
    };

    ZDD<NbInfoFind> nbFindDD(std::max(50 * sphereNumber, nbNumber));
    nbFindDD.clearAll();
    nbFindDD.localID.resize(sphereNumber); // localID and localData means gid and data on local rank
    nbFindDD.localData.resize(sphereNumber);

#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < sphereNumber; i++) {
        nbFindDD.localID[i] = sphere[i].sphereDataIO.gid;

        nbFindDD.localData[i].alpha = sphere[i].sphereDataIO.alpha;
        nbFindDD.localData[i].chebN = sphere[i].myIntegrator->chebN;
        nbFindDD.localData[i].gid = sphere[i].sphereDataIO.gid;
        nbFindDD.localData[i].fdistIndexGlobal = fdistMapRcp->getGlobalElement(sphere[i].fdistIndexLocal);
        // index of assert(nbFindDD.nbLocalData[i].fdistGlobalIndex != Teuchos::OrdinalTraits<int>::invalid());
        nbFindDD.localData[i].direction[0] = sphere[i].sphereDataIO.direction(0);
        nbFindDD.localData[i].direction[1] = sphere[i].sphereDataIO.direction(1);
        nbFindDD.localData[i].direction[2] = sphere[i].sphereDataIO.direction(2);
    }

    nbFindDD.buildIndex();

    // step 2, update nb info for each sphere
    nbFindDD.findID.clear();
    nbFindDD.findData.clear();
    nbFindDD.findID.resize(nbNumber);
    nbFindDD.findData.resize(nbNumber);

#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) { // no openmp
        const int indexBase = nbBufferIndex[i];
        for (int nb = 0; nb < sphere[i].nbBlocks.size(); nb++) {
            // step 1 set a list of nb to get
            nbFindDD.findID[indexBase + nb] = sphere[i].nbBlocks[nb].nbInfo.gid;
        }
    }
    assert(nbFindDD.findID.size() == nbNumber);
    nbFindDD.findData.resize(nbFindDD.findID.size());

    // step 2 load to buffer
    nbFindDD.find();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        // step 3 from buffer to nbBlocks
        const int indexBase = nbBufferIndex[i];
        for (int nb = 0; nb < sphere[i].nbBlocks.size(); nb++) {
            const int nbbufferIndex = indexBase + nb;
            sphere[i].nbBlocks[nb].alpha = nbFindDD.findData[nbbufferIndex].alpha;
            sphere[i].nbBlocks[nb].fdistIndexGlobal = nbFindDD.findData[nbbufferIndex].fdistIndexGlobal;
            sphere[i].nbBlocks[nb].direction[0] = nbFindDD.findData[nbbufferIndex].direction[0];
            sphere[i].nbBlocks[nb].direction[1] = nbFindDD.findData[nbbufferIndex].direction[1];
            sphere[i].nbBlocks[nb].direction[2] = nbFindDD.findData[nbbufferIndex].direction[2];
            const int chebIndex = (nbFindDD.findData[nbbufferIndex].chebN - this->chebIntegrator[0].chebN) /
                                  (this->chebIntegrator[1].chebN - this->chebIntegrator[0].chebN);
            assert(chebIndex >= 0 && chebIndex < chebIntegrator.size());
            sphere[i].nbBlocks[nb].chebIntegrator = &(this->chebIntegrator[chebIndex]);
            assert(sphere[i].nbBlocks[nb].chebIntegrator->chebN == nbFindDD.findData[nbbufferIndex].chebN);
#ifdef ZDDDEBUG
            std::cout << sphere[i].nbBlocks[nb].nbInfo.gid << std::endl;
            std::cout << sphere[i].nbBlocks[nb].nbInfo.pos << std::endl;
            std::cout << sphere[i].nbBlocks[nb].direction << std::endl;
            std::cout << sphere[i].nbBlocks[nb].chebIntegrator->chebN << std::endl;
            std::cout << sphere[i].nbBlocks[nb].fdistIndexGlobal << std::endl;
#endif
        }
    }
    commRcp->barrier();
    std::cout << "finished nb ZDD lookup" << std::endl;

// step 2 update neighbor blocks
#pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].calculateNeighborBlocks(updateNbBlockMat, stkReg);
        // allocate the space for nbRegMat and nbBlock in this routine
    }
    std::cout << "nbBlocks updated" << std::endl;
}

void HydroSphereSystem::locateSphere() {
    // load the spheresendtarget to the findData buffer of sphereGidFindDD
    // pre: valid sphereIO list
    // post: valid sphereGidFindDD.findData list

    sphereGidFindDD.clearAll();
    sphereGidFindDD.localID.resize(sphereIO.size());
    sphereGidFindDD.localData.resize(sphereIO.size());
    for (int i = 0; i < sphereIO.size(); i++) {
        sphereGidFindDD.localID[i] = sphereIO[i].gid;
        sphereGidFindDD.localData[i] = myRank;
    }

    sphereGidFindDD.buildIndex();
    sphereGidFindDD.findID.resize(sphere.size());
    sphereGidFindDD.findData.resize(sphere.size());
    for (int i = 0; i < sphere.size(); i++) {
        sphereGidFindDD.findID[i] = sphere[i].sphereDataIO.gid;
    }

    // load to buffer
    sphereGidFindDD.find();
    // for (auto &srank : sphereGidFindDD.findData) {
    // std::cout << myRank << " " << srank << std::endl;
    // }

    return;
}

void HydroSphereSystem::sortSphere() {
    // sort sphere according to the sphereIO
    // setup an unordered map of gid-index pair

    // pre: valid sphereIO and sphere list
    // post: sphere and sphereIO in the same order

    assert(sphere.size() == sphereIO.size());
    // construct a hash map
    std::unordered_map<int, int> mapGidIndex;
    for (int i = 0; i < sphereIO.size(); i++) {
        mapGidIndex[sphereIO[i].gid] = i;
    }
    // define the compare object
    class objCompare {
      public:
        std::unordered_map<int, int> *mapPtr;
        objCompare(std::unordered_map<int, int> &mapGidIndex) { mapPtr = &mapGidIndex; }
        bool operator()(const RigidSphere &a, const RigidSphere &b) {
            const auto &iteA = mapPtr->find(a.sphereDataIO.gid);
            const auto &iteB = mapPtr->find(b.sphereDataIO.gid);
            if (iteA == mapPtr->end() || iteB == mapPtr->end()) {
                std::cout << "error looking for index\n";
                exit(1);
            }
            const int indexA = iteA->second;
            const int indexB = iteB->second;
            // assert(indexA != indexB);
            return (indexA < indexB);
        }
    };

    // sort
    std::sort(sphere.begin(), sphere.end(), objCompare(mapGidIndex));
    // check
    for (int i = 0; i < sphere.size(); i++) {
        assert(sphere[i].sphereDataIO.gid == sphereIO[i].gid);
    }

    return;
}

void HydroSphereSystem::exchangeSphere() {
    // move spheres according to the location of sphereIO
    // step 1, decide the mpi send target for each sphere
    commRcp->barrier();
    locateSphere();
    assert(sphere.size() == sphereGidFindDD.findData.size());

    // step 2, pack, send, recv
    std::vector<int> nsend(nProcs, 0);
    std::vector<int> nsend_disp(nProcs + 1, 0);
    std::vector<int> nrecv(nProcs, 0);
    std::vector<int> nrecv_disp(nProcs + 1, 0);
    MPI_Request *req_send = new MPI_Request[nProcs];
    MPI_Request *req_recv = new MPI_Request[nProcs];

    // calculate where to go
    int nSphere = sphere.size();
    for (int ip = 0; ip < nSphere; ip++) {
        int srank = sphereGidFindDD.findData[ip];
        if (srank != myRank) {
            nsend[srank]++;
        }
    }
    nsend_disp[0] = 0;
    for (int i = 0; i < nProcs; i++) {
        nsend_disp[i + 1] += nsend_disp[i] + nsend[i];
    }

    // send buffer
    sendBuff.resize(nsend_disp[nProcs]);
    // std::cout << "rank" << myRank << " " << sendBuff.size() << std::endl;
    // *** align send particles on ptcl_send_ *************
    for (int i = 0; i < nProcs; i++) {
        nsend[i] = 0;
    }
    // loop over all particles, leave only those stay in this rank.

    int iloc = 0;
    nSphere = sphere.size();
    assert(sphere.size() == sphereGidFindDD.findData.size());
    // for (int i = 0; i < nSphere; i++) {
    //     std::cout << sphere[i].fdist << std::endl;
    //     std::cout << sphere[i].fdistlast << std::endl;
    // }
    for (int ip = 0; ip < nSphere; ip++) {
        const int srank = sphereGidFindDD.findData[ip];
        if (srank == myRank) {
            // particles[iloc] = particles[ip]; // or use swap
            // std::cout << "swapping rank" << myRank << "iloc" << iloc << std::endl;
            swap(sphere[iloc], sphere[ip]);
            // std::cout << "swapped rank" << myRank << "ip" << ip << std::endl;
            iloc++;
        } else {
            int jloc = nsend[srank] + nsend_disp[srank];
            // std::cout << "rank" << myRank << "srank" << srank << "jloc" << jloc << std::endl;
            // sendBuff[jloc] = particles[ip];
            packRigidSphere(sphere[ip], sendBuff[jloc]);
            // std::cout << "packed data rank" << myRank << "jloc" << jloc << std::endl;
            nsend[srank]++;
        }
    }

    sphere.resize(iloc);
    MPI_Barrier(MPI_COMM_WORLD);
    // receive the number of receive particles
    MPI_Alltoall(nsend.data(), 1, MPI_INT, nrecv.data(), 1, MPI_INT, MPI_COMM_WORLD);
    nrecv_disp[0] = 0;
    for (int i = 0; i < nProcs; i++) {
        nrecv_disp[i + 1] = nrecv_disp[i] + nrecv[i];
    }

    // prepare recv buffer
    recvBuff.resize(nrecv_disp[nProcs]);

    // actual send & recv
    int n_proc_send = 0;
    int n_proc_recv = 0;
    for (int ib = 1; ib < nProcs; ib++) {
        int idsend = (ib + myRank) % nProcs;
        if (nsend[idsend] > 0) {
            int adrsend = nsend_disp[idsend];
            int tagsend = (myRank < idsend) ? myRank : idsend;
            // req_send[n_proc_send++] = MPI::COMM_WORLD.Isend(ptcl_send_.getPointer(adrsend), nsend[idsend],
            // GetDataType<Tptcl>(), idsend, tagsend);
            // send as char type
            int status = MPI_Isend(&(sendBuff[adrsend]), nsend[idsend], GetDataType<SpherePacked>(), idsend, tagsend,
                                   MPI_COMM_WORLD, &(req_send[n_proc_send]));
            if (status != 0) {
                std::cout << "send error" << std::endl;
                exit(1);
            }
            n_proc_send++;
        }
        int idrecv = (nProcs + myRank - ib) % nProcs;
        if (nrecv[idrecv] > 0) {
            int adrrecv = nrecv_disp[idrecv];
            int tagrecv = (myRank < idrecv) ? myRank : idrecv;
            // req_recv[n_proc_recv++] = MPI::COMM_WORLD.Irecv(ptcl_recv_.getPointer(adrrecv), nrecv[idrecv],
            // GetDataType<Tptcl>(), idrecv, tagrecv);
            int status = MPI_Irecv(&(recvBuff[adrrecv]), nrecv[idrecv], GetDataType<SpherePacked>(), idrecv, tagrecv,
                                   MPI_COMM_WORLD, &(req_recv[n_proc_recv]));
            if (status != 0) {
                std::cout << "recv error" << std::endl;
                exit(1);
            }
            n_proc_recv++;
        }
    }
    MPI_Waitall(n_proc_send, req_send, MPI_STATUSES_IGNORE);
    MPI_Waitall(n_proc_recv, req_recv, MPI_STATUSES_IGNORE);

    // step 3, unpack
    // put particles in recvBuffer to the
    // particles.insert(particles.end(), recvBuffer.begin(), recvBuffer.end());
    for (int i = 0; i < recvBuff.size(); i++) {
        sphere.emplace_back();
        unpackRigidSphere(recvBuff[i], sphere.back());
    }
    assert(sphere.size() == sphereIO.size());

    // for (const auto &par : sendBuffer) {
    // std::cout << par.gid << std::endl;
    // }

    // step 4, shift the order to make sure sphere and sphereIO has the same order of gid
    sortSphere();
    assert(sphere.size() == sphereIO.size());
    for (int i = 0; i < sphere.size(); i++) {
        assert(sphere[i].sphereDataIO.gid == sphereIO[i].gid);
    }

    // clean to save memory
    sendBuff.clear();
    recvBuff.clear();
    sphereGidFindDD.clearAll();
    delete[] req_send;
    delete[] req_recv;

    return;
}

void HydroSphereSystem::updatefdistIndex() {
    // pre: valid sphere list
    // post: valid contigMap and valid fdistIndexLocal, valid totalChebN
    int sumChebN = 0;
    const int sphereNumber = sphere.size();
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].fdistIndexLocal = sumChebN;
        sumChebN += 3 * sphere[i].myIntegrator->chebN + 3;
    }

    totalChebN = sumChebN;
    fdistMapRcp =
        Teuchos::rcp(new TMAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), totalChebN, 0, commRcp));
    mobilityMapRcp =
        Teuchos::rcp(new TMAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), sphereNumber * 6, 0, commRcp));

    // initialize temporary calculation space
    fdistRcp = Teuchos::rcp(new TV(fdistMapRcp, false));
    fdistTempSpaceRcp = Teuchos::rcp(new TV(fdistMapRcp, false));

    const int globalSize = fdistMapRcp->getGlobalNumElements();
    std::vector<int> globalIndexOnLocal(globalSize);
#pragma omp parallel for schedule(dynamic, 1024)
    for (int kk = 0; kk < globalSize; kk++) {
        globalIndexOnLocal[kk] = kk;
    }
    fdistFullMapRcp = Teuchos::rcp(new TMAP(globalSize, globalIndexOnLocal.data(), globalSize, 0, commRcp));
    return;
}

void HydroSphereSystem::updateFMMTree() {
    // pre: valid fdistIndexLocal, valid totalChebN, initialized myFMM
    // post: myFMM tree updated, src/trg array with proper size
    src_coord.resize(totalChebN);
    std::cout << "total cheb N" << totalChebN << std::endl;
    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        const int pointsN = sphere[i].myIntegrator->chebN + 1;
        const int fdistIndex = sphere[i].fdistIndexLocal;
        const auto &sphereData = sphere[i].sphereDataIO;
        const auto myIntegrator = sphere[i].myIntegrator;
        for (int j = 0; j < pointsN; j++) {
            src_coord[fdistIndex + 3 * j] =
                (sphereData.center[0] + myIntegrator->points[j] * sphereData.direction[0] * sphereData.alpha);
            src_coord[fdistIndex + 3 * j + 1] =
                (sphereData.center[1] + myIntegrator->points[j] * sphereData.direction[1] * sphereData.alpha);
            src_coord[fdistIndex + 3 * j + 2] =
                (sphereData.center[2] + myIntegrator->points[j] * sphereData.direction[2] * sphereData.alpha);
        }
    }
    trg_coord = src_coord;
    commRcp->barrier();
    myFMM.FMM_UpdateTree(src_coord, trg_coord);
    src_value.resize(totalChebN);
    trg_value.resize(totalChebN);
    std::fill(trg_value.begin(), trg_value.end(), 0.0);
    std::fill(src_value.begin(), src_value.end(), 0.0);
    std::cout << "FMM tree updated" << std::endl;
}

void HydroSphereSystem::dumpSphere(int level) const {
    std::cout << "***********************************" << std::endl;
    std::cout << "dumping sphere" << std::endl;
    std::cout << "rank: " << myRank << " sphere size" << sphere.size() << std::endl;
    for (auto &rigidSphere : sphere) {

        std::cout << "--------------" << std::endl;
        std::cout << "rank: " << myRank << " " << rigidSphere.selfBlock.size() << "\n";
        rigidSphere.sphereDataIO.dumpCout(level);
        if (level > 0) {
            for (int i = 0; i < rigidSphere.nbBlocks.size(); i++) {
                std::cout << "nb gid " << rigidSphere.nbBlocks[i].nbInfo.gid << "nbRegMat size"
                          << rigidSphere.nbBlocks[i].nbRegMat.size() << "nbBlockMat size"
                          << rigidSphere.nbBlocks[i].nbBlockMat.size() << std::endl;
            }
        }
        // std::cout << "vdist" << rigidSphere.vdist << std::endl;
        std::cout << "--------------" << std::endl;
    }
    std::cout << "***********************************" << std::endl;
}

void HydroSphereSystem::dumpSphereIO(int level) const {
    std::cout << "***********************************" << std::endl;
    std::cout << "dumping sphereIO" << std::endl;
    std::cout << "rank " << myRank << " sphereIO size" << sphereIO.size() << std::endl;
    for (auto &sphereIO : sphereIO) {
        std::cout << "--------------" << std::endl;
        std::cout << "rank " << myRank << std::endl;
        sphereIO.dumpCout(level);
        std::cout << "--------------" << std::endl;
    }
    std::cout << "***********************************" << std::endl;
}

void HydroSphereSystem::packRigidSphere(const RigidSphere &rigidSphere, SpherePacked &packedData) {
    // pre: valid two references
    // post:  copy data to the packed struct
    packedData.gid = rigidSphere.sphereDataIO.gid;
    packedData.alpha = rigidSphere.sphereDataIO.alpha; // dimensionless length
    packedData.b = rigidSphere.sphereDataIO.b;         // dimensionless slender parameter
    packedData.meshDelta = rigidSphere.sphereDataIO.b; // largest allowable mesh size
    for (int i = 0; i < 3; i++) {
        // input
        packedData.center[i] = rigidSphere.sphereDataIO.center[i];
        packedData.direction[i] = rigidSphere.sphereDataIO.direction[i];
        packedData.force[i] = rigidSphere.sphereDataIO.force[i];
        packedData.torque[i] = rigidSphere.sphereDataIO.torque[i];
        // output
        packedData.xdot[i] = rigidSphere.sphereDataIO.xdot[i];
        packedData.tdot[i] = rigidSphere.sphereDataIO.tdot[i];
    }

    // the intermediate data
    packedData.chebN = rigidSphere.myIntegrator->chebN;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            packedData.mobMatRot[3 * i + j] = rigidSphere.mobMatRot(i, j);
            packedData.mobMatTrans[3 * i + j] = rigidSphere.mobMatTrans(i, j);
        }
    }
    const int distSize = 3 * rigidSphere.myIntegrator->chebN + 3;
    copyFromEvec(rigidSphere.fdist, &(packedData.fdist[0]), distSize);
    copyFromEvec(rigidSphere.fdistlast, &(packedData.fdistlast[0]), distSize);
    // copyFromEvec(rigidSphere.k0fdist, &(packedData.k0fdist[0]), distSize);
    // copyFromEvec(rigidSphere.Vmdist, &(packedData.Vmdist[0]), distSize);
    // copyFromEvec(rigidSphere.K0Vmdist, &(packedData.K0Vmdist[0]), distSize);
    // copyFromEvec(rigidSphere.vdistself, &(packedData.vdistself[0]), distSize);
    // copyFromEvec(rigidSphere.fdistout, &(packedData.fdistout[0]), distSize);
    // copyFromEvec(rigidSphere.vdist, &(packedData.vdist[0]), distSize);
    // copyFromEvec(rigidSphere.vdistReg, &(packedData.vdistReg[0]), distSize);

    // fill the rest of the array with empty data. unnecessary
    // for (int i = packedData.chebN + 1; i < HYDROMESHMAX * 3; i++) {
    // }

    return;
};

void HydroSphereSystem::unpackRigidSphere(const SpherePacked &packedData, RigidSphere &rigidSphere) {
    // pre: valid two references
    // post: copy packed struct data to rigidSphere, nbBlocks/nbInfos clear, pointer is valid, Evec/Emat with proper
    // size
    rigidSphere.sphereDataIO.gid = packedData.gid;
    rigidSphere.sphereDataIO.alpha = packedData.alpha; // dimensionless length
    rigidSphere.sphereDataIO.b = packedData.b;         // dimensionless slender parameter
    rigidSphere.sphereDataIO.meshDelta = packedData.b; // largest allowable mesh size
    for (int i = 0; i < 3; i++) {
        // input
        rigidSphere.sphereDataIO.center[i] = packedData.center[i];
        rigidSphere.sphereDataIO.direction[i] = packedData.direction[i];
        rigidSphere.sphereDataIO.force[i] = packedData.force[i];
        rigidSphere.sphereDataIO.torque[i] = packedData.torque[i];
        // output
        rigidSphere.sphereDataIO.xdot[i] = packedData.xdot[i];
        rigidSphere.sphereDataIO.tdot[i] = packedData.tdot[i];
    }
    int chebIntIndex = (packedData.chebN / HYDROMESHSTEP) - 1;
    rigidSphere.myIntegrator = &(chebIntegrator[chebIntIndex]);
    assert(rigidSphere.myIntegrator->chebN == packedData.chebN);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rigidSphere.mobMatRot(i, j) = packedData.mobMatRot[3 * i + j];
            rigidSphere.mobMatTrans(i, j) = packedData.mobMatTrans[3 * i + j];
        }
    }

    const int chebN = rigidSphere.myIntegrator->chebN;
    const int distSize = 3 * chebN + 3;
    rigidSphere.fdist.resize(distSize);
    rigidSphere.fdistlast.resize(distSize);
    rigidSphere.allocTemporaryVec(); // set size according to myIntegrator
    copyToEvec(rigidSphere.fdist, &(packedData.fdist[0]), distSize);
    copyToEvec(rigidSphere.fdistlast, &(packedData.fdistlast[0]), distSize);
    // copyToEvec(rigidSphere.k0fdist, &(packedData.k0fdist[0]), distSize);
    // copyToEvec(rigidSphere.Vmdist, &(packedData.Vmdist[0]), distSize);
    // copyToEvec(rigidSphere.K0Vmdist, &(packedData.K0Vmdist[0]), distSize);
    // copyToEvec(rigidSphere.vdistself, &(packedData.vdistself[0]), distSize);
    // copyToEvec(rigidSphere.fdistout, &(packedData.fdistout[0]), distSize);
    // copyToEvec(rigidSphere.vdist, &(packedData.vdist[0]), distSize);
    // copyToEvec(rigidSphere.vdistReg, &(packedData.vdistReg[0]), distSize);

    rigidSphere.selfBlock.resize(distSize, distSize);

    rigidSphere.nbBlocks.clear();
    rigidSphere.nbBlocks.reserve(50);
    // reconstruct the self block and resample mesh

    return;
};

void HydroSphereSystem::saveFdist() {
    // save fdist to fdistlast
    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].fdistlast = sphere[i].fdist;
    }
}

void HydroSphereSystem::reSampleEvecXYZ(Evec &xyzvec, const int newChebN, const ChebNodal *chebInt) { return; }

void HydroSphereSystem::constructDiscretization() {
    // pre: valid sphere list
    // post: valid myIntegrator pointer, valid selfblock, valid dist
    const int sphereNumber = sphere.size();
#pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < sphereNumber; i++) {
        auto &rigidSphere = sphere[i];
        // TODO: resample and allocate space

        printf("rows: %ld\n", rigidSphere.selfBlock.rows());
        rigidSphere.calculateSelfBlock();
        rigidSphere.calculateMobilityMatrix();
    }

    std::cout << "discretization updated" << std::endl;
    return;
}

ChebNodal *HydroSphereSystem::getChebIntegrator(const double &meshDelta, const double &alpha) { return nullptr; }

void HydroSphereSystem::assembleFlowRegMatrix(double droptol) {
    // TODO: assemble close regularization matrix
    return;
}

void HydroSphereSystem::assembleApproxImpSolverMatrix(double droptol) {
    // pre: valid sphere list with nbBlocks.nbBlockMat already calculated
    // post: ApproxImpSolverMatrix assembled, nbBlocks.nbBlockMat resized to zero to save memory access
    // TODO: assemble the close evaluation matrix

    // dumpTCMAT(approxImpSolverMatrixRcp, "approxImpSolverMatrixRcp.mtx");

    return;
}

void HydroSphereSystem::calcSphereMotionWithFdist(SolverType solverType) {
    // pre: valid fdist, vdist etc in sphere list, fdistRcp is synced with this information
    // post: valid xdot/tdot in sphere list
    // Get the number of vectors (columns) in X (and Y).
    // TODO:
}

void HydroSphereSystem::setFdistInitialGuess(Teuchos::RCP<TV> &fdistGuessRcp) {
    // fill fdistGuessRcp with force
    // replace fdistlast with fdist. (save the information for next timestep)
    assert(fdistGuessRcp->getLocalLength() == this->totalChebN);
    // setup initial guess with sphere.fdist, fdistlast
    auto x_2d = fdistGuessRcp->getLocalView<Kokkos::HostSpace>();
    fdistGuessRcp->modify<Kokkos::HostSpace>();

    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        const int fdistIndexLocal = sphere[i].fdistIndexLocal;
        const int chebPoints = sphere[i].myIntegrator->chebN + 1;
        assert(sphere[i].fdist.size() == chebPoints * 3);
        assert(sphere[i].fdist.size() == sphere[i].fdistlast.size());
        for (int j = 0; j < chebPoints * 3; j++) {
            const double x = sphere[i].fdist[j];
            const double y = sphere[i].fdistlast[j];
            x_2d(fdistIndexLocal + j, 0) = 2 * x - y; // linear extrapolation
        }
    }
}
void HydroSphereSystem::solveFullImplicit() {
    // implicitly solve background flow with density
    // pre: valid sphereIO, no repeated gid index. valid input force/torque in sphereIO
    // after: valid velocity/omega in sphereIO

    syncSphereInForceTorque();
    const int sphereNumber = sphere.size();

    // TODO: prepare the linear system
    // initialize b

    // TODO: make an initial guess
    // use data in fdistRcp as initialGuess
    // setting initial is useless with LCP
    // setFdistInitialGuess(fdistSolRcp);

    // #ifdef HYDRODEBUG
    dumpTV(bRcp, "bRcp.mtx");
    dumpTV(fdistSolRcp, "fdistSolRcp.mtx");
    dumpTV(fdistRcp, "fdistRcp.mtx");
    // #endif
    // TODO: setup iterative solver
    if (bRcp->norm2() < 1e-10) {
        // b = 0, no need to solve, return a zero fdist
        fdistSolRcp->putScalar(0);
        fdistRcp->putScalar(0);
        fdistTempSpaceRcp->putScalar(0);
#pragma omp parallel for
        for (int i = 0; i < sphereNumber; i++) {
            sphere[i].clearAll();
        }
        commRcp->barrier();
    } else {
        // assemble hydro operator, the matrix A in Ax=b
        RCP<TOP> AopRcp = rcp(new HydroSphereOperator(this));
        commRcp->barrier();
        // set Belos object
        Belos::SolverFactory<double, TMV, TOP> factory;
        // Make an empty new parameter list.
        Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
        solverParams->set("Num Blocks", 500); // usless for BicgStab. the restart m for GMRES.
        solverParams->set("Maximum Iterations", 5000);
        // solverParams->set("Orthogonalization", "IMGS"); // DGKS, ICGS, IMGS
        solverParams->set("Convergence Tolerance", this->iteres);
        solverParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary);
        // Create the GMRES solver.
        auto solverRCP = factory.create("GMRES", solverParams);
        // BicgSTAB takes less iterations but each iteration takes two Op*x, so no reduction in computation time.
        // the residual is slightly smaller than GMRES
        // left preconditioner with BICGSTAB is not implemented
        auto problemRCP = Teuchos::rcp(new Belos::LinearProblem<TOP::scalar_type, TMV, TOP>(AopRcp, fdistSolRcp, bRcp));
        problemRCP->setProblem(); // necessary

        // get preconditioner
        Teuchos::RCP<TOP> precOp;
        switch (precChoice) {
        case PrecType::ILUT:
            precOp = createILUTPreconditioner(this->approxImpSolverMatrixRcp, 1e-5, 1.0); // droptol,fill
            break;
        case PrecType::RELAXATION:
            precOp = createPlnPreconditioner(this->approxImpSolverMatrixRcp);
            break;
        case PrecType::KINV:
            // dynamic cast from TV to base class TMV
            precOp = createKinvPreconditioner(this->approxImpSolverMatrixRcp, fdistSolRcp);
            break;
        default:
            precOp = createPlnPreconditioner(this->approxImpSolverMatrixRcp);
        }
        problemRCP->setRightPrec(precOp); // or left
        std::cout << solverRCP->description() << std::endl;
        solverRCP->setProblem(problemRCP);
        commRcp->barrier();

        // get solution
        Belos::ReturnType result = solverRCP->solve();
        // apply one more time without the preconditioner, generate the correct vdist, kdist, etc
        AopRcp->apply(*fdistSolRcp, *fdistTempSpaceRcp);
    }

    // TODO: get motion with force
    commRcp->barrier();
    calcSphereMotionWithFdist(SolverType::IMPLICIT);
    // generate the correct, 'unpreconditioned' xdot and tdot in sphere.

    saveFdist();

    // TODO: get result
    syncSphereOut();
    dumpSphere(4);

    return;
}

void HydroSphereSystem::solveExplicitWithBackgroundFlow() {
    // explicitly solve with known background flow

    // step 1, calculate background flow
    // vdist is filled on every sphere without the 1/8pi factor
    syncSphereInForceTorque();

    calcBackgroundFlow(); // one FMM call

    calcSphereMotionWithFdist(SolverType::EXPLICIT);

    updateFdistWithBackgroundflow(); // a small linear system for each sphere ifndef KSETZERO

    saveFdist();

    syncSphereOut();
    return;
}

void HydroSphereSystem::solveMobilityVelocityBackgroundFlow() {
    // with known fdist, calculate the motion when force/torque = 0;
    // WARNING: when calling this for LCP Collision, the sphereIO.force and sphereIO.torque should
    // include all non-collision forces
    syncSphereInForceTorque();
    calcBackgroundFlow();
    calcSphereMotionWithFdist(SolverType::EXPLICIT);

    // allocate a new vector if necessary
    if (mobilityVelocityBackgroundFlowRcp.is_null() == true ||
        mobilityVelocityBackgroundFlowRcp->getMap()->isSameAs(*mobilityMapRcp) != true) {
        mobilityVelocityBackgroundFlowRcp = Teuchos::rcp(new TV(mobilityMapRcp, false));
    }
    assert(mobilityVelocityBackgroundFlowRcp->getMap()->isSameAs(*mobilityMapRcp));

    const int sphereNumber = sphere.size();

    // step  put xdot/tdot into mobilityVelocity
    auto VB_2d = mobilityVelocityBackgroundFlowRcp->getLocalView<Kokkos::HostSpace>();
    mobilityVelocityBackgroundFlowRcp->modify<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        auto &fdata = sphere[i].sphereDataIO;
        // special case for sphere where torque is always perpendicular to direction
        Evec3 omega = fdata.direction.cross(fdata.tdot);
        VB_2d(6 * i + 0, 0) = fdata.xdot[0];
        VB_2d(6 * i + 1, 0) = fdata.xdot[1];
        VB_2d(6 * i + 2, 0) = fdata.xdot[2];
        VB_2d(6 * i + 3, 0) = omega[0];
        VB_2d(6 * i + 4, 0) = omega[1];
        VB_2d(6 * i + 5, 0) = omega[2];
    }

    return;
}

void HydroSphereSystem::updateFdistWithForceTorque() {
    // with force and torque input to sphereIO, calculate the fdist
    // WARNING: when calling this for LCP Collision, the sphereIO.force and sphereIO.torque should
    // include all non-collision forces PLUS collision forces
    syncSphereInForceTorque();

    updateFdistWithBackgroundflow(); // a small linear system for each sphere ifndef KSETZERO
    return;
}

void HydroSphereSystem::updateFdistWithBackgroundflow() {
    // pre valid vdist in sphere list
    // post fdist in sphere list updated

    // TODO: algorithm must be developed
    return;
}

void HydroSphereSystem::applyImpSolverMatrixFdist(const TMV &fdistIn, TMV &fdistOut, double alpha, double beta) {
    // Y := beta*Y + alpha*Op(A)*X
    // the apply function of HydroSphereOperator.
    // pre: all things are ready
    // post: data in sphere list updated, fdistOutRcp updated

    // typedef Teuchos::ScalarTraits<double> STS;
    // Get the number of vectors (columns) in X (and Y).
    // TODO:

    return;
}

void HydroSphereSystem::calcBackgroundFlow() {
    // calculate background flow with info in sphere.fdist
    // pre: valid sphere list, allocated fdistRcpa and fdistTempRcp
    // post: fdistRcp sync with sphere.fdist, sphere.vdist, vdistReg, vdistself updated
    constexpr double pimul8 = (8 * 3.141592653589793238462643383279);

    assert(totalChebN == fdistRcp->getLocalLength());
    assert(totalChebN == fdistTempSpaceRcp->getLocalLength());
    assert(totalChebN == src_value.size());
    assert(totalChebN == src_coord.size());
    assert(totalChebN == trg_value.size());
    assert(totalChebN == trg_coord.size());

    // vdist, vdistself, vdistReg do not include the 1/8pi factor.

    return;
}

void HydroSphereSystem::setPrecType(PrecType choice) { this->precChoice = choice; }

void HydroSphereSystem::calcFlowOnGridWithFdist(double Lscale, double Tscale, double dx, std::string filename,
                                                double xShift, double yShift, double zShift) {
    // pre: valid sphere list and fdist
    // post: pointFlowList hold the flow in the structured grid with predifined max size,
    //      flow dumpped to a file
    // generate the fluid velocity with 3D Cartesian mesh size dx, and write it to filename
    // using Lscale, Tscale as the dimension scale of mesh size.
    // adding a shift if necessary.

    // pay attention that the FMM box size is usually larger than the sphere periodic boundary size
    // except for the TP:PXYZ case. Therefore the shift should be carefully set to match the data

    // step 1 generate cartesian mesh on rank0
    // total mesh size
    const int maxMesh = 50;
    const int NX = std::min((boxHigh[0] - boxLow[0]) / dx, maxMesh * 1.0); // set to 200 costs 8GB memory, too large
    const int NY = std::min((boxHigh[1] - boxLow[1]) / dx, maxMesh * 1.0);
    const int NZ = std::min((boxHigh[2] - boxLow[2]) / dx, maxMesh * 1.0);
    // max 100 points in each dimension
    dx = (boxHigh[0] - boxLow[0]) / (NX);
    const double dy = (boxHigh[1] - boxLow[1]) / (NY);
    const double dz = (boxHigh[2] - boxLow[2]) / (NZ);

    if (myRank == 0) {

        /* TO match VTK point ordering:
         *  # NOTE: VTK expects data in FORTRAN order
         *  The order and number of points must match that specified by the dimensions of the grid.
         *  The point order increases in i fastest (from 0<=i<dims[0]), then j (0<=j<dims[1]), then k (0<=k<dims[2])
         *  where dims[] are the dimensions of the grid in the i-j-k topological directions.
         *  The number of points is dims[0]*dims[1]*dims[2].
         *
         *  The same is true for the cells of the grid.
         *  The order and number of cells must match that specified by the dimensions of the grid.
         *  The cell order increases in i fastest (from 0<=i<(dims[0]-1)), then j (0<=j<(dims[1]-1)),
         *  then k (0<=k<(dims[2]-1)) The number of cells is (dims[0]-1)*(dims[1]-1)*(dims[2]-1).
         * */
        trg_coord.resize(NX * NY * NZ * 3);
        trg_value.resize(NX * NY * NZ * 3);
#pragma omp parallel for
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    // const Evec3 pos = Evec3(i * dx, j * dy, k * dz) + boxLow;
                    trg_coord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 0] = (i + 0.5) * dx + boxLow[0];
                    trg_coord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 1] = (j + 0.5) * dy + boxLow[1];
                    trg_coord[(i + j * (NX) + k * (NX) * (NY)) * 3 + 2] = (k + 0.5) * dz + boxLow[2];
                }
            }
        }
        std::fill(trg_value.begin(), trg_value.end(), 0.0);
    } else {
        trg_coord.clear();
        trg_value.clear();
    }
    pointFlowList.clear();
    commRcp->barrier();

    // step 2 set up FMM tree and calc velocity, on every node
    const int sphereN = sphere.size();
    // update totalChebN and fdistIndexLocal
    int sumChebN = 0;
    const int sphereNumber = sphere.size();
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].fdistIndexLocal = sumChebN;
        sumChebN += 3 * sphere[i].myIntegrator->chebN + 3;
    }
    this->totalChebN = sumChebN;
    src_coord.resize(this->totalChebN);

#pragma omp parallel for
    for (int i = 0; i < sphereN; i++) {
        const int pointsN = sphere[i].myIntegrator->chebN + 1;
        const int fdistIndex = sphere[i].fdistIndexLocal;
        const auto &fib = sphere[i];
        const auto &fdata = sphere[i].sphereDataIO;
        for (int j = 0; j < pointsN; j++) {
            src_coord[fdistIndex + 3 * j + 0] =
                (fdata.center[0] + fib.myIntegrator->points[j] * fdata.direction[0] * fdata.alpha);
            src_coord[fdistIndex + 3 * j + 1] =
                (fdata.center[1] + fib.myIntegrator->points[j] * fdata.direction[1] * fdata.alpha);
            src_coord[fdistIndex + 3 * j + 2] =
                (fdata.center[2] + fib.myIntegrator->points[j] * fdata.direction[2] * fdata.alpha);
        }
    }

    commRcp->barrier();
    myFMM.FMM_UpdateTree(src_coord, trg_coord);

    // setup source values
    src_value.resize(this->totalChebN);
#pragma omp parallel for
    for (int i = 0; i < sphereN; i++) {
        const int pointsN = sphere[i].myIntegrator->chebN + 1;
        const int fdistIndex = sphere[i].fdistIndexLocal;
        const auto &fib = sphere[i];
        const auto &fdata = sphere[i].sphereDataIO;
        for (int j = 0; j < pointsN; j++) {
            src_value[fdistIndex + 3 * j + 0] = sphere[i].fdist[3 * j + 0] * fib.myIntegrator->weights[j] * fdata.alpha;
            src_value[fdistIndex + 3 * j + 1] = sphere[i].fdist[3 * j + 1] * fib.myIntegrator->weights[j] * fdata.alpha;
            src_value[fdistIndex + 3 * j + 2] = sphere[i].fdist[3 * j + 2] * fib.myIntegrator->weights[j] * fdata.alpha;
        }
    }
    commRcp->barrier();
    myFMM.FMM_TreeClear();
    myFMM.FMM_Evaluate(trg_value, trg_value.size() / 3, &src_value);

    if (myRank == 0) {
        // step 3 save data
        pointFlowList.resize(NX * NY * NZ);
        const int pointFlowNumber = pointFlowList.size();
#pragma omp parallel for
        for (int i = 0; i < pointFlowNumber; i++) {
            // convert from dimensionless to dimensional
            // the FMM call results contain the 1/8pi factor.
            pointFlowList[i].pos[0] = trg_coord[3 * i + 0] * Lscale;
            pointFlowList[i].pos[1] = trg_coord[3 * i + 1] * Lscale;
            pointFlowList[i].pos[2] = trg_coord[3 * i + 2] * Lscale;
            pointFlowList[i].vel[0] = trg_value[3 * i + 0] * Lscale / Tscale;
            pointFlowList[i].vel[1] = trg_value[3 * i + 1] * Lscale / Tscale;
            pointFlowList[i].vel[2] = trg_value[3 * i + 2] * Lscale / Tscale;
        }

        // step 4 dump
        std::string vtk(".vtk");
        filename = filename + vtk;
        FILE *fdump = fopen(filename.c_str(), "w");
        // VTK file header
        fprintf(fdump, "# vtk DataFile Version 3.0\n");
        fprintf(fdump, "flow velocity\n");
        fprintf(fdump, "ASCII\n");
        fprintf(fdump, "DATASET STRUCTURED_POINTS\n");
        fprintf(fdump, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
        fprintf(fdump, "ORIGIN %f %f %f\n", (boxLow[0] + dx * 0.5) * Lscale, (boxLow[1] + dy * 0.5) * Lscale,
                (boxLow[2] + dz * 0.5) * Lscale);
        fprintf(fdump, "SPACING %f %f %f\n", dx * Lscale, dy * Lscale, dz * Lscale);

        // VTK point properties
        fprintf(fdump, "POINT_DATA %d\n", NX * NY * NZ);
        fprintf(fdump, "SCALARS Velocity float 3\n"); //	SCALARS dataName dataType numComp
        fprintf(fdump, "LOOKUP_TABLE DEFAULT\n");     // (NOT optional) default look up table

        for (int i = 0; i < pointFlowNumber; i++) {
            fprintf(fdump, "%.6e %.6e %.6e\n", pointFlowList[i].vel[0], pointFlowList[i].vel[1],
                    pointFlowList[i].vel[2]);
        }

        fclose(fdump);
    }

    pointFlowList.clear();

    commRcp->barrier();
    trg_value.clear();
    trg_coord.clear();
    src_value.clear();
    src_coord.clear();
    myFMM.FMM_TreeClear();
}

Teuchos::RCP<TOP> &HydroSphereSystem::getImplicitMobilityOperator(const double Lscale, const double Tscale,
                                                                  const double Escale) {
    // get the mobility matrix
    mobilityOperatorRcp = Teuchos::rcp(new MobilityOperator(this, Lscale, Tscale, Escale));
    return mobilityOperatorRcp;
}

Teuchos::RCP<TV> &HydroSphereSystem::getMobilityVelocityBackgroundFlow(const double Lscale, const double Tscale,
                                                                       const double Escale) {
    const double Uscale = Lscale / Tscale;
    const double Oscale = 1 / Tscale;

    // step 3 put xdot/tdot into mobilityVelocity
    auto VB_2d = mobilityVelocityBackgroundFlowRcp->getLocalView<Kokkos::HostSpace>();
    mobilityVelocityBackgroundFlowRcp->modify<Kokkos::HostSpace>();

    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        // special case for sphere where torque is always perpendicular to direction
        VB_2d(6 * i + 0, 0) = VB_2d(6 * i + 0, 0) * Uscale; // velocity
        VB_2d(6 * i + 1, 0) = VB_2d(6 * i + 1, 0) * Uscale;
        VB_2d(6 * i + 2, 0) = VB_2d(6 * i + 2, 0) * Uscale;
        VB_2d(6 * i + 3, 0) = VB_2d(6 * i + 3, 0) * Oscale; // omega
        VB_2d(6 * i + 4, 0) = VB_2d(6 * i + 4, 0) * Oscale;
        VB_2d(6 * i + 5, 0) = VB_2d(6 * i + 5, 0) * Oscale;
    }

#ifdef HYDRODEBUG
    std::cout << "mobilityVelocityBackgroundFlow constructed" << mobilityVelocityBackgroundFlowRcp->description()
              << std::endl;
    dumpTV(mobilityVelocityBackgroundFlowRcp, "mobilityVelocityBackgroundFlow.mtx");
#endif

    return this->mobilityVelocityBackgroundFlowRcp;
}

Teuchos::RCP<const TMAP> &HydroSphereSystem::getFdistMap() { return this->fdistMapRcp; }

Teuchos::RCP<const TMAP> &HydroSphereSystem::getMobilityMap() { return this->mobilityMapRcp; }
}
