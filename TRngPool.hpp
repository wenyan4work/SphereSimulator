
#ifndef TRNGPOOL_H_
#define TRNGPOOL_H_

#include "mpi.h"
#include "omp.h"

#include "trng/config.hpp"
#include "trng/lcg64_shift.hpp"
#include "trng/lognormal_dist.hpp"
#include "trng/mrg5.hpp"
#include "trng/normal_dist.hpp"
#include "trng/uniform01_dist.hpp"

#include <iostream>
#include <memory>

class TRngPool {
  private:
    int myRank;
    int nProcs;
    int nThreads;

    // typedef trng::lcg64_shift myEngineType;
    typedef trng::mrg5 myEngineType;
    std::vector<std::unique_ptr<myEngineType>> rngEngineThreadsPtr;

    trng::uniform01_dist<double> u01;
    trng::normal_dist<double> n01;

    std::unique_ptr<trng::lognormal_dist<double>> lnDistPtr;

  public:
    void setLogNormalParameters(double mu, double sigma) {
        lnDistPtr.reset(new trng::lognormal_dist<double>(mu, sigma));
    }

    explicit TRngPool(int seed = 0) : n01(0, 1) {

        myRank = 0;
        nProcs = 1;
        nThreads = 1;
        int mpiInitFlag;
        MPI_Initialized(&mpiInitFlag);
        if (mpiInitFlag > 0) {
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
        } else {
            myRank = 0;
            nProcs = 1;
        }
        nThreads = omp_get_max_threads();
        myEngineType rngEngine;
        rngEngine.seed(static_cast<unsigned long>(seed));

        if (nProcs > 1) {
            rngEngine.split(nProcs, myRank);
        }

        rngEngineThreadsPtr.resize(nThreads);
#pragma omp parallel for num_threads(nThreads)
        for (int i = 0; i < nThreads; i++) {
            // a copy of engine for each thread
            rngEngineThreadsPtr[i].reset(new myEngineType());
            *rngEngineThreadsPtr[i] = rngEngine;
            // split
            rngEngineThreadsPtr[i]->split(nThreads, i);
        }
        setLogNormalParameters(1.0, 1.0);
    };

    ~TRngPool() = default;

    TRngPool(const TRngPool &) = delete;
    TRngPool &operator=(const TRngPool &) = delete;

    inline double getU01(int threadId) { return u01(*rngEngineThreadsPtr[threadId]); }
    inline double getN01(int threadId) { return n01(*rngEngineThreadsPtr[threadId]); }
    inline double getLN(int threadId) { return (*lnDistPtr)(*rngEngineThreadsPtr[threadId]); }

    inline double getU01() {
        const int threadId = omp_get_thread_num();
        return getU01(threadId);
    }

    inline double getN01() {
        const int threadId = omp_get_thread_num();
        return getN01(threadId);
    }

    inline double getLN() {
        const int threadId = omp_get_thread_num();
        return getLN(threadId);
    }
};

#endif