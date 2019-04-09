/*
    Created by Dhairya Malhotra
    Modified by Wen Yan
    2017-2018
*/

#ifndef _NEAR_INTERAC_HPP_
#define _NEAR_INTERAC_HPP_

#include "sctl/sctl.hpp"

// Interface of OBJ:
// const double * Coord() const;
// const double Rad() const;
// void Pack(std::vector<char> &buff) const ;
// void Unpack(const std::vector<char> &buff) ;
// optional:  void CopyFromFull(const FullType &) ;

template <class Real, int DIM>
struct Pair {
    using Long = sctl::Long;

    Long trg_idx, src_idx; // the ID in the local trgNear and srcNear containers, not the global ID
    Real srcShift[DIM];    // the shift added to src.coord() where the pair is detected

    int operator<(const Pair<Real, DIM> &p1) const {
        return std::pair<Long, Long>(trg_idx, src_idx) < std::pair<Long, Long>(p1.trg_idx, p1.src_idx);
    }
};

template <class Real, int DIM>
class NearInteraction {
    using MID = sctl::Morton<DIM>;
    using Long = sctl::Long;

    // internal data for each object, trivially copyable
    struct ObjData {
        int operator<(const ObjData &p1) const { return mid < p1.mid; }
        MID mid;   // Morton ID
        Long Rglb; // global sorted unique id

        Real rad;
        Real coord[DIM]; // coord supplied, not necessarily in the original box
    };

  public:
    NearInteraction() { Init(); }

    NearInteraction(sctl::Comm comm) : comm_(comm) { Init(); }

    void SetPeriodLength(sctl::Integer d, Real len) {
        SCTL_ASSERT(d < DIM);
        period_length[d] = len;
    }

    template <class SrcObj, class TrgObj>
    void SetupRepartition(const std::vector<SrcObj> &src_vec, const std::vector<TrgObj> &trg_vec);

    template <class SrcObj, class TrgObj>
    void SetupNearInterac(const std::vector<SrcObj> &src_vec, const std::vector<TrgObj> &trg_vec);

    const std::vector<Pair<Real, DIM>> &GetInteractionList() const { return trg_src_pair; }

    template <class ObjType>
    void ForwardScatterSrc(const std::vector<ObjType> &in, std::vector<ObjType> &out) const;

    template <class ObjType>
    void ForwardScatterTrg(const std::vector<ObjType> &in, std::vector<ObjType> &out) const;

    template <class ObjType>
    void ReverseScatterTrg(const std::vector<ObjType> &in, std::vector<ObjType> &out) const;

    void Barrier() { comm_.Barrier(); }

  private:
    void Init() {
        // default is non-periodic
        for (sctl::Integer i = 0; i < DIM; i++) {
            period_length[i] = 0;
            period_length0[i] = 0;
        }
    }

    template <class ObjType>
    void ForwardScatter(const std::vector<ObjType> &in_vec, std::vector<ObjType> &out_vec,
                        const sctl::Vector<Long> &recv_idx) const;

    template <class ObjType>
    void ReverseScatter(const std::vector<ObjType> &in_vec, std::vector<ObjType> &out_vec,
                        const sctl::Vector<Long> &send_idx) const;

    sctl::Comm comm_;

    /***********************************
     * These are set in SetupPartition()
     ***********************************/
    sctl::Integer depth;  // the depth is set such that a src-trg pair must appear within two neighboring leaf octants

    sctl::Vector<Long> TRglb, SRglb;      // globally ordered, locally sorted sequential internal ID
    sctl::Vector<ObjData> SData_, TData_; // Globally sorted ObjData

    Real period_length[DIM];
    Real period_length0[DIM]; // peroidic length and scaled periodic length
    // Real s = sctl::pow<Real>(2, depth);
    // period_length0[i] = std::floor((period_length[i] / BBox.L) * s) / s;

    /***********************************
     * These are set in SetupNearInterac()
     ***********************************/
    sctl::Vector<Pair<Real, DIM>> TSPair;
    std::vector<Pair<Real, DIM>> trg_src_pair;
};

#include "NearInteraction.txx"

#endif //_NEAR_INTERAC_HPP_
