#ifndef INTERACTIONMANAGER_HPP
#define INTERACTIONMANAGER_HPP

#include <cstdio>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "NearInteraction.hpp"

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

  public:
    // constructor
    InteractionManager(std::vector<SrcObjType> *srcVecPtr_, std::vector<TrgObjType> *trgVecPtr_) {
        srcVecPtr = srcVecPtr_;
        trgVecPtr = trgVecPtr_;
        if (srcVecPtr == trgVecPtr) {
            printf("Working for src == trg\n");
        }
    }

    // default copy
    InteractionManager(const InteractionManager &) = default;
    InteractionManager(InteractionManager &&) = default;
    InteractionManager &operator=(const InteractionManager &) = default;
    InteractionManager &operator=(InteractionManager &&) = default;

    ~InteractionManager(){
        // do not delete the two pointers
    };

    // get a new NearInteraction object to be used for partition or interaction
    std::shared_ptr<NearInteraction<Real, DIM>> getNewNearInteraction() {
        sctl::Comm sctlcomm;
#ifdef SCTL_HAVE_MPI
        sctlcomm = MPI_COMM_WORLD;
#endif
        return std::make_shared<NearInteraction<Real, DIM>>(sctlcomm);
    }

    // load balancing depending on near interaction. this is supposed to be not frequently called.
    // Note the special case TrgObjType=SrcObjType, and trgVecPtr=srcVecPtr
    void partitionObject(std::shared_ptr<NearInteraction<Real, DIM>> &nearInteracPtr) {
        // after partition, the containers pointed by srcVecPtr and trgVecPtr are modified

        // Repartition data
        // Setup for repartition
        nearInteracPtr->template SetupRepartition<SrcObjType, TrgObjType>(*srcVecPtr, *trgVecPtr);
        // nearInteracPtr->Barrier();
        MPI_Barrier(MPI_COMM_WORLD);

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
                             std::vector<SrcEssType> &srcEssVec, std::vector<TrgEssType> &trgEssVec) const {

        // Compute near interactions
        // Setup for near interaction
        nearInteracPtr->SetupNearInterac(srcEssVec, trgEssVec);
        // nearInteracPtr->Barrier();
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    template <class SrcEssType, class TrgEssType, class InteractionFunctorType>
    void calcNearInteraction(const std::shared_ptr<NearInteraction<Real, DIM>> &nearInteracPtr,
                             std::vector<SrcEssType> &srcEssVec, std::vector<TrgEssType> &trgEssVec,
                             InteractionFunctorType &interactor) const {
        // nearInteracPtr->Barrier();
        MPI_Barrier(MPI_COMM_WORLD);

        // It is the user's duty to make sure that the nearInteracPtr, srcEssVec, and trgEssVec
        // here is equal to the nearInteracPtr used in the setupNear() funciton

        // Following code can repeat multiple times without calling Setup again as long as
        // the position and shape of the sources and the targets do not change.

        // Forward scatter
        std::vector<SrcEssType> srcNear;
        std::vector<TrgEssType> trgNear;
        nearInteracPtr->template ForwardScatterSrc<SrcEssType>(srcEssVec, srcNear);
        nearInteracPtr->template ForwardScatterTrg<TrgEssType>(trgEssVec, trgNear);

        // trg_src_interac is sorted according to trg id
        // This is the internal id maintained by the object pointed by nearInteracPtr
        // This is different from the user defined data field gid in Full type or Ess type

        const auto &trg_src_interac = nearInteracPtr->GetInteractionList(); // Get interaction list

        // Debug
#ifdef DEBUGINTERACT
        std::cout << "-----------" << std::endl;
        for (auto &v : trg_src_interac) {
            std::cout << v.first << " " << v.second << std::endl;
        }
        std::cout << "-----------" << std::endl;
#endif

        long N = trg_src_interac.size();

        if (N == 0) {
            // no pairs to process
            return;
        }

        std::vector<std::pair<int, int>> trgIdIndex;
        trgIdIndex.reserve(trg_src_interac.size());
        trgIdIndex.emplace_back(std::pair<int, size_t>(trg_src_interac[0].first, 0));
        for (size_t i = 1; i < N; i++) {
            if (trg_src_interac[i].first != trg_src_interac[i - 1].first) {
                // the index within N of a new trg ID.
                trgIdIndex.emplace_back(std::pair<int, size_t>(trg_src_interac[i].first, i));
            }
        }
        trgIdIndex.emplace_back(std::pair<int, size_t>(static_cast<int>(UINTMAX_MAX), N)); // a tail to simplify index

        const int ntrg = trgIdIndex.size() - 1; // the last element is only for index bound
#pragma omp parallel for
        for (size_t k = 0; k < ntrg; k++) {
            // start index
            size_t lb = trgIdIndex[k].second;
            // end index
            size_t ub = trgIdIndex[k + 1].second;

// real work
#pragma omp simd
            for (size_t i = lb; i < ub; i++) {
                TrgEssType &t = trgNear[trg_src_interac[i].first];
                SrcEssType &s = srcNear[trg_src_interac[i].second];
                // ( compute interaction between t and s )
                // potentially define this as #pragma omp declare simd
                interactor(t, s);
            }
        }

#ifdef DEBUGINTERACT
        std::cout << "show pair" << std::endl;
        for (auto &pair : trg_src_interac) {
            std::cout << trgNear[pair.first].gid << " " << srcNear[pair.second].gid << std::endl;
        }

        std::cout << "show index" << std::endl;
        for (auto &p : trgIdIndex) {

            std::cout << p.first << " " << p.second << std::endl;
        }
#endif

        // Reverse scatter
        // nearInteracPtr->Barrier();
        MPI_Barrier(MPI_COMM_WORLD);
        nearInteracPtr->template ReverseScatterTrg<TrgEssType>(trgNear, trgEssVec);

        return;
    }

  private:
};

#endif