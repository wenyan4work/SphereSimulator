#ifndef _NEAR_INTERAC_HPP_
#define _NEAR_INTERAC_HPP_

#include "Util/sctl.hpp"

template <class Real, int DIM> class NearInteraction {
    typedef sctl::Morton<DIM> MID;
    typedef sctl::Long Long;
    struct ObjData {
      int operator<(const ObjData& p1) const { return mid < p1.mid; }
      MID mid;
      Long Rglb;

      Real rad;
      Real coord[DIM];
      //sctl::StaticArray<Real,DIM> coord; // not trivially copyable
    };

  public:
    NearInteraction() {}

    NearInteraction(sctl::Comm comm) : comm_(comm) {}

    template <class SrcObj, class TrgObj> void SetupRepartition(const std::vector<SrcObj>& src_vec, const std::vector<TrgObj>& trg_vec);

    template <class SrcObj, class TrgObj> void SetupNearInterac(const std::vector<SrcObj>& src_vec, const std::vector<TrgObj>& trg_vec);

    const std::vector<std::pair<Long,Long>>& GetInteractionList() const { return trg_src_pair; }

    template <class ObjType> void ForwardScatterSrc(const std::vector<ObjType>& in, std::vector<ObjType>& out) const;

    template <class ObjType> void ForwardScatterTrg(const std::vector<ObjType>& in, std::vector<ObjType>& out) const;

    template <class ObjType> void ReverseScatterTrg(const std::vector<ObjType>& in, std::vector<ObjType>& out) const;

    void Barrier(){comm_.Barrier();}

  private:
    template <class ObjType> void ForwardScatter(const std::vector<ObjType>& in_vec, std::vector<ObjType>& out_vec, const sctl::Vector<Long>& recv_idx) const;

    template <class ObjType> void ReverseScatter(const std::vector<ObjType>& in_vec, std::vector<ObjType>& out_vec, const sctl::Vector<Long>& send_idx) const;

    sctl::Comm comm_;
    sctl::Vector<Long> TRglb, SRglb;
    sctl::Vector<std::pair<Long,Long>> TSPair;
    std::vector<std::pair<Long, Long>> trg_src_pair;

    Real period_length;
    sctl::Integer depth;
    sctl::Vector<ObjData> SData_, TData_; // Globally sorted
};

#include "NearInteraction.txx"

#endif //_NEAR_INTERAC_HPP_
