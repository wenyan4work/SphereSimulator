#ifndef INTERACTIONMANAGER_HPP
#define INTERACTIONMANAGER_HPP

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

#include "NearInteraction.hpp"

#include <mpi.h>
#include <omp.h>

// 1. given object, repartition
// 2. given essential type, interaction functor, compute near interaction

// for 1 species:
// SrcObjType == TrgObjType
// trgVecPtr = srcVecPtr
// SrcEssType can be equal or not to TrgEsstype

// for 2 species:
// SrcObjType != TrgObjType
// trgVecPtr != srcVecPtr
// SrcEssType can be equal or not to TrgEsstype

template <class Real, int DIM, class SrcObjType, class TrgObjType>
class InteractionManager {
    std::vector<SrcObjType> *srcVecPtr;
    std::vector<TrgObjType> *trgVecPtr;

    std::array<Real, DIM> boxLow;
    std::array<Real, DIM> boxHigh;
    std::array<bool, DIM> pbcFlag; // 111 for periodic in xyz, etc

    // WARNING: create or copy this too many times will cause fatal error:
    //   Too many communicators (0/16384 free on this process; ignore_id=0)
    sctl::Comm sctlcomm;

  public:
    // constructor
    InteractionManager(std::vector<SrcObjType> *srcVecPtr_, std::vector<TrgObjType> *trgVecPtr_) {
        srcVecPtr = srcVecPtr_;
        trgVecPtr = trgVecPtr_;
        if (srcVecPtr == trgVecPtr) {
            printf("Working for src == trg\n");
        }
#ifdef SCTL_HAVE_MPI
        sctlcomm = MPI_COMM_WORLD;
#endif

        for (int i = 0; i < DIM; i++) {
            boxLow[i] = 0;
            boxHigh[i] = 0;
            pbcFlag[i] = 0;
        }
    }

    // default copy
    InteractionManager(const InteractionManager &) = default;
    InteractionManager(InteractionManager &&) = default;
    InteractionManager &operator=(const InteractionManager &) = default;
    InteractionManager &operator=(InteractionManager &&) = default;

    ~InteractionManager() = default; // do not delete the two pointers

    // set PBC
    // If PBC is activated, it is the user's job to guarantee that the coordinates are within the box set by PBC
    // There is no coordinate check of that
    void setPBCBox(std::array<bool, DIM> &pbcFlag_, std::array<Real, DIM> &boxLow_, std::array<Real, DIM> &boxHigh_) {
        pbcFlag = pbcFlag_;
        for (int i = 0; i < DIM; i++) {
            if (boxLow_[i] >= boxHigh_[i]) {
                printf("box low must < box high \n");
                exit(1);
            }
            boxLow[i] = boxLow_[i];
            boxHigh[i] = boxHigh_[i];
        }
    }

    // a helper routine for PBC
    void fitInPeriodicBound(Real &x, const Real &lb, const Real &ub) {
        // put the periodic image of x in [lb,ub)
        const Real L = ub - lb;
        while (x >= ub) {
            x -= L;
        }
        while (x < lb) {
            x += L;
        }
    }

    std::shared_ptr<NearInteraction<Real, DIM>> getNewNearInteraction() {
        // get a new NearInteraction object to be used for partition or interaction
        auto p = std::make_shared<NearInteraction<Real, DIM>>(sctlcomm);
        for (int i = 0; i < DIM; i++) {
            if (pbcFlag[i]) {
                p->SetPeriodLength(i, boxHigh[i] - boxLow[i]);
            } else {
                p->SetPeriodLength(i, 0);
            }
        }
        return p;
    }

    // load balancing depending on near interaction. this is supposed to be not frequently called.
    // Note the special case TrgObjType=SrcObjType, and trgVecPtr=srcVecPtr
    void partitionObject(std::shared_ptr<NearInteraction<Real, DIM>> &nearInteracPtr) {
        // after partition, the containers pointed by srcVecPtr and trgVecPtr are modified

        // Repartition data
        // Setup for repartition
        nearInteracPtr->template SetupRepartition<SrcObjType, TrgObjType>(*srcVecPtr, *trgVecPtr);
        MPI_Barrier(sctlcomm.GetMPI_Comm());

        // Distribute source and target vectors
        std::vector<SrcObjType> srcNew;
        nearInteracPtr->template ForwardScatterSrc<SrcObjType>(*srcVecPtr, srcNew);
        srcVecPtr->swap(srcNew);

        if (trgVecPtr != srcVecPtr) {
            std::vector<TrgObjType> trgNew;
            nearInteracPtr->template ForwardScatterTrg<TrgObjType>(*trgVecPtr, trgNew);
            trgVecPtr->swap(trgNew);
        }
    }

    template <class SrcEssType, class TrgEssType>
    void setupEssVec(std::vector<SrcEssType> &srcEssVec, std::vector<TrgEssType> &trgEssVec) const {
        // copy from full type to ess type
        const int nsrc = srcVecPtr->size();
        srcEssVec.resize(nsrc);
#pragma omp parallel for
        for (int i = 0; i < nsrc; i++) {
            srcEssVec[i].CopyFromFull((*srcVecPtr)[i]);
        }

        const int ntrg = trgVecPtr->size();
        trgEssVec.resize(ntrg);
#pragma omp parallel for
        for (int i = 0; i < ntrg; i++) {
            trgEssVec[i].CopyFromFull((*trgVecPtr)[i]);
        }
    }

    template <class SrcEssType, class TrgEssType>
    void setupNearInteractor(const std::shared_ptr<NearInteraction<Real, DIM>> &nearInteracPtr,
                             std::vector<SrcEssType> &srcEssVec, std::vector<TrgEssType> &trgEssVec) {

        // Compute near interactions
        // Setup for near interaction

        nearInteracPtr->SetupNearInterac(srcEssVec, trgEssVec);
        MPI_Barrier(sctlcomm.GetMPI_Comm());
        return;
    }

    template <class SrcEssType, class TrgEssType, class InteractionFunctorType>
    void calcNearInteraction(const std::shared_ptr<NearInteraction<Real, DIM>> &nearInteracPtr,
                             std::vector<SrcEssType> &srcEssVec, std::vector<TrgEssType> &trgEssVec,
                             InteractionFunctorType &interactor) {
        MPI_Barrier(sctlcomm.GetMPI_Comm());

        // It is the user's duty to make sure that the nearInteracPtr, srcEssVec, and trgEssVec
        // here is equal to the nearInteracPtr used in the setupNear() funciton

        // Following code can repeat multiple times without calling Setup again as long as
        // the position and shape of the sources and the targets do not change.

        // Forward scatter
        // printf("scatter start\n");
        std::vector<SrcEssType> srcNear;
        std::vector<TrgEssType> trgNear;
        nearInteracPtr->template ForwardScatterSrc<SrcEssType>(srcEssVec, srcNear);
        nearInteracPtr->template ForwardScatterTrg<TrgEssType>(trgEssVec, trgNear);
        // printf("scatter complete\n");

        // trg_src_interac is sorted according to trg id
        // This is the internal id maintained by the object pointed by nearInteracPtr
        // This is different from the user defined data field gid in Full type or Ess type

        const auto &trg_src_interac = nearInteracPtr->GetInteractionList(); // Get interaction list

        const long N = trg_src_interac.size();

        if (N == 0) {
            // no pairs to process
            // do nothing
        } else {
            // printf("calc work division start\n");
            std::vector<std::pair<size_t, size_t>> trgIdIndex;
            trgIdIndex.reserve(trgNear.size());
            trgIdIndex.emplace_back(std::pair<int, size_t>(trg_src_interac[0].trg_idx, 0));
            for (size_t i = 1; i < N; i++) {
                if (trg_src_interac[i].trg_idx != trg_src_interac[i - 1].trg_idx) {
                    // the index within N of a new trg ID.
                    trgIdIndex.emplace_back(std::pair<size_t, size_t>(trg_src_interac[i].trg_idx, i));
                }
            }

#ifdef DEBUGINTERACT
            std::cout << "---------------" << std::endl;
            std::cout << "show pair" << std::endl;
            for (auto &pair : trg_src_interac) {
                std::cout << pair.trg_idx << " " << pair.src_idx << std::endl;
            }

            std::cout << "show index" << std::endl;
            for (auto &p : trgIdIndex) {

                std::cout << p.first << " " << p.second << std::endl;
            }
            std::cout << "---------------" << std::endl;
#endif

            // printf("calc work division complete\n");
            const int ntrg = trgIdIndex.size(); // the last element is only for index bound
#pragma omp parallel for
            for (size_t k = 0; k < ntrg; k++) {
                size_t lb = trgIdIndex[k].second;                                             // start index
                size_t ub = k < ntrg - 1 ? trgIdIndex[k + 1].second : trg_src_interac.size(); // end index
                for (size_t i = lb; i < ub; i++) {
                    // real work
                    const auto &trg_idx = trg_src_interac[i].trg_idx;
                    const auto &src_idx = trg_src_interac[i].src_idx;
                    SCTL_ASSERT(trg_idx < trgNear.size());
                    SCTL_ASSERT(src_idx < srcNear.size());
                    TrgEssType &t = trgNear[trg_idx];
                    SrcEssType &s = srcNear[src_idx];
                    std::array<Real, DIM> srcShift;
                    for (int k = 0; k < DIM; k++) {
                        srcShift[k] = trg_src_interac[i].srcShift[k];
                    }
                    // ( compute interaction between t and s )
                    // potentially define this as #pragma omp declare simd
#ifdef DEBUGINTERACT
                    std::cout << trg_idx << " " << src_idx << std::endl;
#endif
                    interactor(t, s, srcShift);
                }
            }
        }
        // std::cout << "pair interactor complete" << std::endl;

        // Reverse scatter
        MPI_Barrier(sctlcomm.GetMPI_Comm());
        nearInteracPtr->template ReverseScatterTrg<TrgEssType>(trgNear, trgEssVec);

        return;
    }
};

#endif