#ifndef SYNCMANAGER_HPP
#define SYNCMANAGER_HPP

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <vector>

#include <mpi.h>

#include "Buffer.hpp"

// helper class to broadcast and reduce replicated objects through all process
// required interface of T:
// T.    void Pack(std::vector<char> &buff) const
// T.    void Unpack(const std::vector<char> &buff)
// T.    void Reduce(const T &)
template <class T>
class SyncManager {

    std::vector<T> *objVecPtr;

  public:
    SyncManager(std::vector<T> *objVecPtr_) : objVecPtr(objVecPtr_), rank(0), nProcs(1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    };
    ~SyncManager(){};

    // pack, broadcast from rank 0, unpack
    // precondition: *objVecPtr on ranks other than 0 are cleared
    // postcondition: *objVecPtr on every rank synced with that on rank 0
    void sync() {
        if (nProcs == 1) {
            return;
        }
        int error = 0; // mpi error code
        // step 1 pack al objects.
        std::vector<T> &objVec = *objVecPtr;

        std::vector<int> objSize;
        std::vector<char> objBuffer;

        if (rank == 0) {
            packAll(objVec, objSize, objBuffer);
        } else {
            objVec.clear();
        }

        int nobj = objVec.size();
        assert(objSize.size() == nobj);

        // step 2, broadcast objSize, objBuf
        int size[2] = {objSize.size(), objBuffer.size()};
        // int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
        error = MPI_Bcast(&(size[0]), 2, MPI_INT, 0, MPI_COMM_WORLD);
        assert(error == 0);

        if (rank != 0) {
            objSize.resize(size[0]);
            objBuffer.resize(size[1]);
        }
        error = MPI_Bcast(objSize.data(), size[0], MPI_INT, 0, MPI_COMM_WORLD);
        assert(error == 0);
        error = MPI_Bcast(objBuffer.data(), size[1], MPI_CHAR, 0, MPI_COMM_WORLD);
        assert(error == 0);

        // step 3, unpack and fill objVec on all ranks other than 0
        // this function returns earlier on rank 0 than other ranks
        if (rank != 0) {
            unpackAll(objVec, objSize, objBuffer);
        }
    };

    // pack, copy to rank 9, unpack, reduce every object
    // precondition: *objVecPtr (and each obj in it) on every rank has the same size as rank 0,
    // postcondition: * objVecPtr transferred to rank 0, reduced on rank 0
    // not using the ring-reducing due to complicated pack and unpack overhead
    // objVecPtr must contain the same number and order of objects.
    // same object on different rank may have different size
    void reduce() {
        if (nProcs == 1) {
            return;
        }
        int error = 0; // mpi error code
        std::vector<int> objSize;
        std::vector<char> objBuffer;
        std::vector<int> objSizeAll;
        std::vector<char> objBufferAll;
        // step 1, pack on every rank != 0, allocate buffer on rank=0
        if (rank != 0) {
            packAll(*objVecPtr, objSize, objBuffer);
        } else {
            // objSize.size() is the same on every rank, = objVecPtr->size() on rank 0
            objSize.resize(objVecPtr->size(), 0); // no packed data for rank 0. just occupy space
            objBuffer.resize(0);
            objSizeAll.resize(objVecPtr->size() * nProcs);
        }
        // printf("packed\n");
        assert(objSize.size() == objVecPtr->size());
        // gather size and prepare buffer on rank 0
        // int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
        //  MPI_Datatype recvtype, int root, MPI_Comm comm)
        // int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
        // const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm)

        error = MPI_Gather(objSize.data(), objSize.size(), MPI_INT, objSizeAll.data(), objSize.size(), MPI_INT, 0,
                           MPI_COMM_WORLD);
        // printf("gather\n");
        assert(error == 0);
        // allocate buffer at rank 0
        std::vector<int> recvcounts;
        std::vector<int> displs;
        if (rank == 0) {
            const int nobj = objVecPtr->size();
            assert(nobj == objSize.size());
            assert(nobj * nProcs == objSizeAll.size());
            int bufSize = std::accumulate(objSizeAll.cbegin(), objSizeAll.cend(), 0);
            // printf("buffsize %d\n", bufSize);
            objBufferAll.resize(bufSize);

            recvcounts.resize(nProcs, 0);
            displs.resize(nProcs, 0);
            for (int i = 0; i < nProcs; i++) {
                recvcounts[i] =
                    std::accumulate(objSizeAll.cbegin() + i * nobj, objSizeAll.cbegin() + i * nobj + nobj, 0);
            }
            // for (auto &i : objSizeAll) {
            //     printf("objSizeAll %d\n", i);
            // }
            // for (auto &i : recvcounts) {
            //     printf("recv %d\n", i);
            // }
            std::partial_sum(recvcounts.cbegin(), recvcounts.cend() - 1, displs.begin() + 1);
            // for (auto &i : displs) {
            //     printf("displ %d\n", i);
            // }
        }
        error = MPI_Gatherv(objBuffer.data(), objBuffer.size(), MPI_CHAR, objBufferAll.data(), recvcounts.data(),
                            displs.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
        // printf("gatherv\n");
        assert(error == 0);

        // step 2, unpack and reduce on rank 0
        if (rank == 0) {
            const int nobj = objVecPtr->size();
            // // the location where the buf for obj i from rank j begins in objBufferAll
            std::vector<int> objBuffIndex(objSizeAll.size(), 0);
            std::partial_sum(objSizeAll.cbegin(), objSizeAll.cend() - 1, objBuffIndex.begin() + 1);
            // for (auto &i : objBuffIndex) {
            //     printf("objBuffIndex %d\n", i);
            // }
            // printf("objBufferAll size %d\n", objBufferAll.size());

#pragma omp parallel for
            for (int i = 0; i < nobj; i++) {
                auto &target = (*objVecPtr)[i];
                for (int j = 1; j < nProcs; j++) {
                    // unpack the object i from rank j and construct an temporary object
                    std::vector<char> buf;
                    buf.clear();
                    int bufSize = objSizeAll[nobj * j + i];
                    int bufIndex = objBuffIndex[nobj * j + i];
                    // printf("%d,%d\n", bufSize, bufIndex);
                    std::copy_n(objBufferAll.cbegin() + bufIndex, bufSize, std::back_inserter(buf));
                    // reduce
                    T source;
                    source.Unpack(buf);
                    target.Reduce(source);
                }
            }
        }
    };

  private:
    int rank = 0;
    int nProcs = 1;

    void packAll(const std::vector<T> &objVec, std::vector<int> &objSize, std::vector<char> &objBuffer) {
        // pack from objVec to objSize and objBuffer
        std::vector<char> buf;
        const int nobj = objVec.size();
        for (int i = 0; i < nobj; i++) {
            buf.clear();
            auto &obj = objVec[i];
            obj.Pack(buf);
            int nchar = buf.size();
            objSize.push_back(nchar);
            std::copy(buf.cbegin(), buf.cend(), std::back_inserter(objBuffer));
        }
    }

    void unpackAll(std::vector<T> &objVec, const std::vector<int> &objSize, const std::vector<char> &objBuffer) {
        const int nobj = objSize.size();
        objVec.resize(nobj);
        std::vector<int> objBufferIndex(objSize.size(), 0);
        std::partial_sum(objSize.cbegin(), objSize.cend(), objBufferIndex.begin() + 1);
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < nobj; i++) {
            // get the index of the ith object
            const int index = objBufferIndex[i];
            const int size = objSize[i];
            // a local buffer
            std::vector<char> buf;
            std::copy(objBuffer.cbegin() + index, objBuffer.cbegin() + index + size, std::back_inserter(buf));
            // unpack
            objVec[i].Unpack(buf);
        }
    }
};

#endif