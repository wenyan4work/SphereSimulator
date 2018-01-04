// std
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <string>
#include <utility>

// helper routines
#include "Preconditioner.hpp"
#include "ZDD.hpp"

// namespace HydroSphere
#include "SPHSphereSystem.hpp"

// for debug. NDEBUG defined in pvfmm.hpp
#undef NDEBUG
#include <cassert>

// max N, periodic type, fmm parameter, iterative solver parameter
SPHSphereSystem::SPHSphereSystem(const Evec3 &boxLow_, const Evec3 &boxHigh_, const int &approxNumberSphere,
                                 const int &multiOrder, const FMM_Wrapper::PAXIS &pbc)
    : myFMM(multiOrder, 1000, 0, pbc), commRcp(Tpetra::DefaultPlatform::getDefaultPlatform().getComm()),
      myRank(commRcp->getRank()), nProcs(commRcp->getSize()), sphereGidFindDD(1000), boxLow(boxLow_),
      boxHigh(boxHigh_) {
    // initialization order determined by order in declaration

    //	Kokkos::initialize();

    // initialize FMM object
    myFMM.FMM_SetBox(boxLow[0], boxHigh[0], boxLow[1], boxHigh[1], boxLow[2], boxHigh[2]);

    // Reserve space
    srcValue.reserve(approxNumberSphere * multiOrder * multiOrder * 2);
    trgValue.reserve(approxNumberSphere * multiOrder * multiOrder * 2);
    srcCoord.reserve(approxNumberSphere * multiOrder * multiOrder * 2);
    trgCoord.reserve(approxNumberSphere * multiOrder * multiOrder * 2);

    commRcp->barrier();
    printf("SPHSphereSystem initialized\n");
}

void SPHSphereSystem::addSphere(const SphereIO &newSphere) {
    // precondition: the gid of newSphere must be unique, and valid (>0)
    // postcondition: sphere self data is valid, sphere neighbor blocks are not valid
    sphereIO.emplace_back(newSphere);
    sphere.emplace_back(newSphere);
    return;
}

void SPHSphereSystem::getReadyForOperator() {
    commRcp->barrier();
    // normalize orientation
    const int sphereNumber = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < sphereNumber; i++) {
        sphere[i].orientation.normalize();
    }

    exchangeSphere(); // locateSphere() called in this
    updateNeighborSphere();
    // updateMap();

    return;
}

void SPHSphereSystem::locateSphere() {
    // determine the mpi send target for each sphere with gid i
    // source: current rank of sphere object with gid i
    // destination: the rank of sphereIO object with gid i
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
        sphereGidFindDD.findID[i] = sphere[i].gid;
    }

    // load to buffer
    sphereGidFindDD.find();
    // for (auto &srank : sphereGidFindDD.findData) {
    // std::cout << myRank << " " << srank << std::endl;
    // }

    return;
}

void SPHSphereSystem::sortSphere() {
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
        bool operator()(const Sphere &a, const Sphere &b) {
            const auto &iteA = mapPtr->find(a.gid);
            const auto &iteB = mapPtr->find(b.gid);
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
        assert(sphere[i].gid == sphereIO[i].gid);
    }

    return;
}

void SPHSphereSystem::exchangeSphere() {
    bool debug = false;
#ifndef DEDEBUG
    debug = true;
#endif
    // move spheres according to the location of sphereIO
    // step 1, decide the mpi send target for each sphere
    commRcp->barrier();
    locateSphere();
    const auto &sendDest = sphereGidFindDD.findData;
    if (sphere.size() != sendDest.size()) {
        printf("locate sphere error\n");
    }

    // step 2, pack, send, recv
    // number of spheres
    std::vector<int> nSphSend(nProcs, 0);
    std::vector<int> nSphSendDisp(nProcs + 1, 0);
    std::vector<int> nSphRecv(nProcs, 0);
    std::vector<int> nSphRecvDisp(nProcs + 1, 0);
    // number of bytes
    std::vector<int> nByteSend(nProcs, 0);
    std::vector<int> nByteSendDisp(nProcs + 1, 0);
    std::vector<int> nByteRecv(nProcs, 0);
    std::vector<int> nByteRecvDisp(nProcs + 1, 0);
    // send buffer for each rank
    std::vector<Buffer> sendBuffer(nProcs);
    // each rank has its own send buffer.
    // most would be empty later (no sphere to be sent)

    MPI_Request *req_send = new MPI_Request[nProcs];
    MPI_Request *req_recv = new MPI_Request[nProcs];

    // calculate where to go
    int nSphere = sphere.size();
    for (int ip = 0; ip < nSphere; ip++) {
        int sendRank = sendDest[ip];
        if (sendRank != myRank) {
            nSphSend[sendRank]++;
        }
    }
    nSphSendDisp[0] = 0;
    std::partial_sum(nSphSend.cbegin(), nSphSend.cend(), nSphSendDisp.begin() + 1);

    // fill send buffer
    // loop over all particles, leave only those stay in this rank.
    nSphere = sphere.size();
    assert(sphere.size() == sendDest.size());
    int iloc = 0;
    for (int ip = 0; ip < nSphere; ip++) {
        const int sendRank = sendDest[ip];
        if (sendRank == myRank) {
            // ip should be reserved on this rank.
            // swap ip to the position of iloc, then iloc self increment
            swap(sphere[iloc], sphere[ip]);
            iloc++;
        } else {
            // ip should be sent to
            // int jloc = nsend[sendRank] + nsend_disp[sendRank];
            // std::cout << "rank" << myRank << "srank" << srank << "jloc" << jloc << std::endl;
            // sendBuff[jloc] = particles[ip];
            sphere[ip].pack(sendBuffer[sendRank]);
        }
    }

    sphere.resize(iloc); // discard all spheres to be sent
    // find size in Byte
    for (int i = 0; i < nProcs; i++) {
        nByteSend[i] = sendBuffer[i].getSize();
    }
    nByteSendDisp[0] = 0;
    std::partial_sum(nByteSend.cbegin(), nByteSend.cend(), nByteSendDisp.begin() + 1);

    // output for debug
    if (debug == true) {
        for (int i = 0; i < nProcs; i++) {
            printf("send rank %d, send spheres %d, bytes %d\n", myRank, nSphSend[i], nByteSend[i]);
            printf("send rank %d, disp spheres %d, bytes %d\n", myRank, nSphSendDisp[i], nByteSendDisp[i]);
        }
        printf("send rank %d, disp spheres %d, bytes %d\n", myRank, nSphSendDisp[nProcs], nByteSendDisp[nProcs]);
    }

    commRcp->barrier();

    // receive the number of recv spheres
    MPI_Alltoall(nSphSend.data(), 1, MPI_INT, nSphRecv.data(), 1, MPI_INT, MPI_COMM_WORLD);
    nSphRecvDisp[0] = 0;
    std::partial_sum(nSphRecv.cbegin(), nSphRecv.cend(), nSphRecvDisp.begin() + 1);
    // receive the number of total bytes
    MPI_Alltoall(nByteSend.data(), 1, MPI_INT, nByteRecv.data(), 1, MPI_INT, MPI_COMM_WORLD);
    nByteRecvDisp[0] = 0;
    std::partial_sum(nByteRecv.cbegin(), nByteRecv.cend(), nByteRecvDisp.begin() + 1);

    // output for debug
    if (debug == true) {
        for (int i = 0; i < nProcs; i++) {
            printf("recv rank %d, recv spheres %d, bytes %d\n", myRank, nSphRecv[i], nByteRecv[i]);
            printf("recv rank %d, disp spheres %d, bytes %d\n", myRank, nSphRecvDisp[i], nByteRecvDisp[i]);
        }
        printf("recv rank %d, disp spheres %d, bytes %d\n", myRank, nSphRecvDisp[nProcs], nByteRecvDisp[nProcs]);
    }

    // prepare recv buffer
    std::vector<char> recvBufferRaw(nByteRecvDisp[nProcs], 0);

    // actual send & recv, in asynchronized mode
    int n_proc_send = 0;
    int n_proc_recv = 0;

    for (int ib = 1; ib < nProcs; ib++) {
        // every rank send to myRank+1, recv from myRank-1
        // tag is set to the smaller of the send&recv ranks

        int sendRank = (ib + myRank) % nProcs;
        if (sendBuffer[sendRank].getSize() > 0) {
            int tagsend = (myRank < sendRank) ? myRank : sendRank;
            // send as char type
            int status = MPI_Isend(sendBuffer[sendRank].getPtr(), sendBuffer[sendRank].getSize(), MPI_CHAR, sendRank,
                                   tagsend, MPI_COMM_WORLD, &(req_send[n_proc_send]));
            if (status != 0) {
                printf("send error %d status from rank %d to rank %\n", status, myRank, sendRank);
                exit(1);
            }
            n_proc_send++;
        }
        int recvRank = (nProcs + myRank - ib) % nProcs;
        if (nByteRecv[recvRank] > 0) {
            int adrrecv = nByteRecvDisp[recvRank];
            int tagrecv = (myRank < recvRank) ? myRank : recvRank;
            // recv as char type
            int status = MPI_Irecv(&(recvBufferRaw[adrrecv]), nByteRecv[recvRank], MPI_CHAR, recvRank, tagrecv,
                                   MPI_COMM_WORLD, &(req_recv[n_proc_recv]));
            if (status != 0) {
                printf("recv error %d status for rank %d from rank %\n", status, myRank, recvRank);
                exit(1);
            }
            n_proc_recv++;
        }
    }
    MPI_Waitall(n_proc_send, req_send, MPI_STATUSES_IGNORE);
    MPI_Waitall(n_proc_recv, req_recv, MPI_STATUSES_IGNORE);

    // step 3, unpack
    // put particles in recvBuffer to the
    Buffer recvBuffer(recvBufferRaw);
    while (recvBuffer.getReadPos() < recvBuffer.getSize()) {
        sphere.emplace_back();
        sphere.back.unpack(recvBuffer);
    }
    if (sphere.size() != sphereIO.size()) {
        printf("recv size error\n");
    }
    if (debug) {
        dumpSphere();
        dumpSphereIO();
    }

    // step 4, shift the order to make sure sphere and sphereIO has the same order of gid
    sortSphere();
    assert(sphere.size() == sphereIO.size());
    for (int i = 0; i < sphere.size(); i++) {
        if (sphere[i].gid != sphereIO[i].gid) {
            printf("error in pos %d\n", i);
        }
    }

    // clean to save memory
    sphereGidFindDD.clearAll();
    delete[] req_send;
    delete[] req_recv;

    return;
}

void SPHSphereSystem::updateNeighborSphere() {
    // pre: valid sphereDataIO.nbInfos, valid fdistIndex, valid fdistMapRcp
    // post: valid sphere.nbBlocks
    const int sphereNumber = sphere.size();
    std::vector<int> nbBufferIndex(sphereNumber + 1); // last = total nb number
    nbBufferIndex[0] = 0;
    for (int i = 1; i < sphereNumber + 1; i++) {
        nbBufferIndex[i] = nbBufferIndex[i - 1] + sphereIO[i - 1].sphNeighborIO.size();
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
    std::cout << "nbBlocks updated" << std::endl;
}

void SPHSphereSystem::dumpSphere() const {
    printf("***********************************\n");
    printf("dumping sphere\n");
    printf("rank %d, Sphere size %d\n", myRank, sphere.size());
    for (const auto &sp : sphere) {
        printf("--------------\n");
        printf("rank %d\n", myRank);
        sp.dumpSphere();
        sp.dumpNeighbor();
        for (const auto &l : sp.sphLayer) {
            sp.dumpLayer(l.first);
        }
        printf("--------------\n");
    }
    std::cout << "***********************************" << std::endl;
}

void SPHSphereSystem::dumpSphereIO() const {
    printf("***********************************\n");
    printf("dumping sphereIO\n");
    printf("rank %d, SphereIO size %d\n", myRank, sphereIO.size());
    for (const auto &spIO : sphereIO) {
        printf("--------------\n");
        printf("rank %d\n", myRank);
        spIO.dumpSphere();
        printf("--------------\n");
    }
    printf("***********************************\n");
}
