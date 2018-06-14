#include <omp.h>
template <class Real, int DIM> template <class SrcObj, class TrgObj> void NearInteraction<Real, DIM>::SetupRepartition(const std::vector<SrcObj>& src_vec, const std::vector<TrgObj>& trg_vec) {
  sctl::StaticArray<Long, 2> Nloc, Rglb;
  { // Set Nloc, Rglb
    Nloc[0] = src_vec.size();
    Nloc[1] = trg_vec.size();
    comm_.Scan<Long>(Nloc, Rglb, 2, sctl::Comm::CommOp::SUM);
    Rglb[0] -= Nloc[0];
    Rglb[1] -= Nloc[1];
  }

  sctl::Vector<Real> SRad(Nloc[0]), TRad(Nloc[1]);
  sctl::Vector<Real> SCoord(Nloc[0] * DIM), TCoord(Nloc[1] * DIM);
  #pragma omp parallel for
  for (sctl::Long i = 0; i < Nloc[0]; i++) { // Set SCoord, SRad
    for (sctl::Integer k = 0; k < DIM; k++) {
      SCoord[i * DIM + k] = src_vec[i].Coord()[k];
    }
    SRad[i] = src_vec[i].Rad();
  }
  #pragma omp parallel for
  for (sctl::Long i = 0; i < Nloc[1]; i++) { // Set TCoord, TRad
    for (sctl::Integer k = 0; k < DIM; k++) {
      TCoord[i * DIM + k] = trg_vec[i].Coord()[k];
    }
    TRad[i] = trg_vec[i].Rad();
  }

  //---------------------------------------------------------------------

  struct {
    Real L;
    sctl::StaticArray<Real, DIM> X;
    // Bounding box : [X, X+L)
  } BBox;
  { // Determine bounding box: BBox.X, BBox.L
    sctl::StaticArray<Real, DIM> X0;
    sctl::StaticArray<Real, DIM> X1;
    { // determine local bounding box
      { // Init X0, X1 <-- Broadcast one source or target point to all processors
        sctl::Integer src_proc, src_proc_ = (SCoord.Dim() || TCoord.Dim() ? comm_.Rank() : 0);
        comm_.Allreduce<sctl::Integer>(sctl::Ptr2ConstItr<sctl::Integer>(&src_proc_,1), sctl::Ptr2Itr<sctl::Integer>(&src_proc,1), 1, sctl::Comm::CommOp::MAX);

        if (comm_.Rank() == src_proc) {
          sctl::Iterator<Real> X_ = SCoord.begin();
          if (TCoord.Dim()) X_ = TCoord.begin();
          SCTL_ASSERT(SCoord.Dim() || TCoord.Dim());
          comm_.Allreduce<Real>(X_, X0, DIM, sctl::Comm::CommOp::SUM);
        } else {
          sctl::StaticArray<Real, DIM> X_;
          for (sctl::Integer k = 0; k < DIM; k++) X_[k] = 0;
          comm_.Allreduce<Real>(X_, X0, DIM, sctl::Comm::CommOp::SUM);
        }
        for (sctl::Integer k = 0; k < DIM; k++) X1[k] = X0[k];
      }
      auto determine_bbox = [&](const sctl::Vector<Real>& coord, sctl::Iterator<Real> X0, sctl::Iterator<Real> X1) {
        Long N = coord.Dim() / DIM;
        if(DIM==2){
          Real X00=X0[0],X01=X0[1],X10=X1[0],X11=X1[1];
          #pragma omp parallel for reduction(min:X00,X01) reduction(max:X10,X11)
          for (sctl::Long i = 0; i < N; i++) {
            X00 = std::min(X00, coord[i * DIM + 0]);
            X01 = std::min(X01, coord[i * DIM + 1]);
            X10 = std::max(X10, coord[i * DIM + 0]);
            X11 = std::max(X11, coord[i * DIM + 1]);
          }
          X0[0]=X00;
          X0[1]=X01;
          X1[0]=X10;
          X1[1]=X11;
        }else if(DIM==3){
          Real X00=X0[0],X01=X0[1],X02=X0[2],X10=X1[0],X11=X1[1],X12=X1[2];
          #pragma omp parallel for reduction(min:X00,X01,X02) reduction(max:X10,X11,X12)
          for (sctl::Long i = 0; i < N; i++) {
            X00 = std::min(X00, coord[i * DIM + 0]);
            X01 = std::min(X01, coord[i * DIM + 1]);
            X02 = std::min(X02, coord[i * DIM + 2]);
            X10 = std::max(X10, coord[i * DIM + 0]);
            X11 = std::max(X11, coord[i * DIM + 1]);
            X12 = std::max(X12, coord[i * DIM + 2]);
          }
          X0[0]=X00;
          X0[1]=X01;
          X0[2]=X02;
          X1[0]=X10;
          X1[1]=X11;
          X1[2]=X12;
        }else{
          for (sctl::Long i = 0; i < N; i++) {
            for (sctl::Integer k = 0; k < DIM; k++) {
              X0[k] = std::min(X0[k], coord[i * DIM + k]);
              X1[k] = std::max(X1[k], coord[i * DIM + k]);
            }
          }
        }

      };
      determine_bbox(SCoord, X0, X1);
      determine_bbox(TCoord, X0, X1);
    }

    BBox.L = 0;
    auto& X0glb = BBox.X;
    sctl::StaticArray<Real, DIM> X1glb;
    { // determine global bounding box
      comm_.Allreduce<Real>(X0, X0glb, DIM, sctl::Comm::CommOp::MIN);
      comm_.Allreduce<Real>(X1, X1glb, DIM, sctl::Comm::CommOp::MAX);
      for (sctl::Integer k = 0; k < DIM; k++) { // Set BBox.L = bbox length
        BBox.L = std::max(BBox.L, X1glb[k] - X0glb[k]);
      }
      BBox.L *= 1.01; // so that points are in [0,L) instead of [0,L]
    }
  }
  { // determine tree depth
    depth = 0;
    Real max_rad_glb, max_rad = 0.0,  oct_size = 1.0;
    #pragma omp parallel for reduction(max:max_rad)
    for (sctl::Long i = 0; i < SRad.Dim(); i++) max_rad = std::max(max_rad, SRad[i]);
    #pragma omp parallel for reduction(max:max_rad)
    for (sctl::Long i = 0; i < TRad.Dim(); i++) max_rad = std::max(max_rad, TRad[i]);
    comm_.Allreduce(sctl::Ptr2ConstItr<Real>(&max_rad,1), sctl::Ptr2Itr<Real>(&max_rad_glb,1), 1, sctl::Comm::CommOp::MAX);
    while (0.5 * oct_size * BBox.L > 2 * max_rad_glb) {
      oct_size *= 0.5;
      depth++;
    }
  }
  { // Set period_length0 <-- period_length / BBox.L
    Real s = sctl::pow<Real>(2, depth);
    for (sctl::Integer i = 0; i < DIM; i++) {
      period_length0[i] = std::floor((period_length[i] / BBox.L) * s) / s;
      if (period_length[i] > 0) SCTL_ASSERT(period_length0[i] > 0);
    }
  }

  sctl::Vector<ObjData> SData(Nloc[0]), TData(Nloc[1]);
  { // Set SData, TData
    sctl::StaticArray<Real,DIM> s;
    for (sctl::Integer d = 0; d < DIM; d++) {
      s[d] = 1.0 / BBox.L;
      if (period_length[d] > 0) {
        s[d] = period_length0[d] / period_length[d];
      }
    }

    #pragma omp parallel for
    for (sctl::Long i = 0; i < SData.Dim(); i++) {
    sctl::StaticArray<Real, DIM> X;
      SData[i].rad = SRad[i];
      for (sctl::Integer k = 0; k < DIM; k++) {
        Real XX = SCoord[i * DIM + k];
        SData[i].coord[k] = XX;
        if (period_length[k] > 0) { // Map XX to unit box
          sctl::Long q = std::floor(XX / period_length[k]);
          XX = (XX - q * period_length[k]) * s[k];
        } else {
          XX = (XX - BBox.X[k]) * s[k];
        }
        X[k] = XX;
      }
      SData[i].mid = MID((sctl::ConstIterator<Real>)X, depth);
      SData[i].Rglb = Rglb[0] + i;
    }
    #pragma omp parallel for
    for (sctl::Long i = 0; i < TData.Dim(); i++) {
    sctl::StaticArray<Real, DIM> X;
      TData[i].rad = TRad[i];
      for (sctl::Integer k = 0; k < DIM; k++) {
        Real XX = TCoord[i * DIM + k];
        TData[i].coord[k] = XX;
        if (period_length[k] > 0) { // Map XX to unit box
          sctl::Long q = std::floor(XX / period_length[k]);
          XX = (XX - q * period_length[k]) * s[k];
        } else {
          XX = (XX - BBox.X[k]) * s[k];
        }
        X[k] = XX;
      }
      TData[i].mid = MID((sctl::ConstIterator<Real>)X, depth);
      TData[i].Rglb = Rglb[1] + i;
    }
  }

  { // Set SData_, TData_
    comm_.HyperQuickSort(TData, TData_);
    comm_.HyperQuickSort(SData, SData_);
    SCTL_ASSERT(TData_.Dim());
    auto m0 = TData_[0];
    comm_.PartitionS(SData_, m0);
    comm_.PartitionS(TData_, m0);
  }

  { // Set TRglb, SRglb
    TRglb.ReInit(TData_.Dim());
    SRglb.ReInit(SData_.Dim());
    #pragma omp parallel for
    for (Long i = 0; i < TData_.Dim(); i++) TRglb[i] = TData_[i].Rglb;
    #pragma omp parallel for
    for (Long i = 0; i < SData_.Dim(); i++) SRglb[i] = SData_[i].Rglb;
  }
  sctl::omp_par::merge_sort(TRglb.begin(), TRglb.end());
  sctl::omp_par::merge_sort(SRglb.begin(), SRglb.end());
}

template <class Real, sctl::Integer DIM> static void NbrList(sctl::Vector<sctl::Morton<DIM>>& nbrlst, const sctl::Morton<DIM>& m0, sctl::Integer depth, sctl::ConstIterator<Real> period_length) {
  typedef typename sctl::Morton<DIM>::UINT_T UINT_T;
  const UINT_T N0 = (((UINT_T)1) << depth);
  Real s = ((Real)1) / N0;

  sctl::StaticArray<Real, (1<<((2*DIM)))*DIM> c;
  m0.Coord((sctl::Iterator<Real>)c);
  sctl::Integer len = 1;
  for (sctl::Integer d = 0; d < DIM; d++) { // Set c
    Real c0 = c[d] + 0.5 * s;
    sctl::Integer new_len = 0;
    for (sctl::Integer k = 0; k < 3; k++) {
      Real c0_new = c0 + (k-1) * s;
      if (period_length[d] > 0) {
        sctl::Long q = std::floor(c0_new / period_length[d]);
        c0_new = (c0_new - q * period_length[d]);
      }
      if (c0_new >= 0 && c0_new < 1) {
        for (sctl::Integer i = 0; i < len * DIM; i++) c[(new_len    ) * DIM + i] = c[i];
        for (sctl::Integer i = 0; i < len      ; i++) c[(new_len + i) * DIM + d] = c0_new;
        new_len += len;
      }
    }
    len = new_len;
  }

  nbrlst.ReInit(len);
  for (sctl::Integer i = 0; i < len; i++) { // Set nbrlst
    nbrlst[i] = sctl::Morton<DIM>(c + i * DIM, depth);
  }
}

template <class Real, int DIM> template <class SrcObj, class TrgObj> void NearInteraction<Real, DIM>::SetupNearInterac(const std::vector<SrcObj>& src_vec, const std::vector<TrgObj>& trg_vec) {
  SetupRepartition<SrcObj, TrgObj>(src_vec, trg_vec);
  // printf("repartition setup in nearinterac\n");

  sctl::Vector<ObjData> SData__; // With ghost sources
  { // sctl::Communicate ghost source data
    auto& Data = SData_;
    Long N = Data.Dim();

    sctl::Integer rank = comm_.Rank();
    sctl::Integer np = comm_.Size();
    sctl::Vector<MID> mins(np);
    { // Set mins
      MID m = MID().Next();
      if (!rank) m = MID().DFD(depth);
      if (TData_.Dim()) m = std::min(m, TData_[0].mid);
      if (SData_.Dim()) m = std::min(m, SData_[0].mid);
      comm_.Allgather(sctl::Ptr2ConstItr<MID>(&m,1), 1, mins.begin(), 1);
    }

    sctl::Vector<Long> usr_cnt, usr_dsp, usr_proc;
    { // Determine source user-count and user-processes
      usr_cnt.ReInit(N);
      usr_dsp.ReInit(N);
      usr_cnt = 0;

      std::vector<sctl::Vector<Long>> pusr_tmp_vec(N);
      #pragma omp parallel for
      for (sctl::Long i = 0; i < N; i++) {
        auto & pusr_tmp=pusr_tmp_vec[i];
        pusr_tmp.ReInit(0);
        sctl::Vector<MID> nbrlst;
        NbrList(nbrlst, Data[i].mid, depth, (sctl::ConstIterator<Real>)period_length0);
        for (sctl::Integer k = 0; k < nbrlst.Dim(); k++) {
          Long puser = std::lower_bound(mins.begin(), mins.end(), nbrlst[k].DFD()) - mins.begin() - 1;
          SCTL_ASSERT(0<=puser);
          SCTL_ASSERT(puser<np);
          if (puser != rank) {
            bool insert_puser = true;
            for (sctl::Integer j = 0; j < pusr_tmp.Dim(); j++) {
              if (pusr_tmp[j] == puser) {
                insert_puser = false;
                break;
              }
            }
            if (insert_puser) {
              pusr_tmp.PushBack(puser);
            }
          }
        }
        // usr_dsp[i] = usr_proc.Dim();
        usr_cnt[i] += pusr_tmp.Dim();
      }
      usr_dsp[0]=0;
      for(int i=1;i<N;i++){
        usr_dsp[i]=usr_cnt[i-1]+usr_dsp[i-1];
      }
      for(auto &pusr_tmp:pusr_tmp_vec)
        for (auto &p:pusr_tmp)
          usr_proc.PushBack(p);
    }

    sctl::Vector<std::pair<Long,ObjData>> usr_data_pair;
    if (N) {
      usr_data_pair.ReInit(usr_dsp[N - 1] + usr_cnt[N - 1]);
      #pragma omp parallel for
      for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Long j = 0; j < usr_cnt[i]; j++) {
          Long idx = usr_dsp[i] + j;
          usr_data_pair[idx].first = usr_proc[idx];
          usr_data_pair[idx].second = Data[i];
        }
      }
      sctl::omp_par::merge_sort(usr_data_pair.begin(), usr_data_pair.end());
    }

    { // Set SData__
      sctl::Vector<Long> scnt(np), sdsp(np);
      scnt = 0;
      sdsp = 0;
      sctl::Vector<ObjData> sbuff(usr_data_pair.Dim());
      #pragma omp parallel for
      for (sctl::Long i = 0; i < usr_data_pair.Dim(); i++) {
        Long p = usr_data_pair[i].first;
        sbuff[i] = usr_data_pair[i].second;
        #pragma omp atomic update
        scnt[p]++;
      }
      sctl::omp_par::scan(scnt.begin(), sdsp.begin(), np);

      sctl::Vector<Long> rcnt(np), rdsp(np);
      rcnt = 0;
      rdsp = 0;
      comm_.Alltoall(scnt.begin(), 1, rcnt.begin(), 1);
      SCTL_ASSERT(rcnt[rank] == 0);
      rcnt[rank] = N;
      sctl::omp_par::scan(rcnt.begin(), rdsp.begin(), np);
      SData__.ReInit(rdsp[np - 1] + rcnt[np - 1]);
      rcnt[rank] = 0;

      comm_.Alltoallv(sbuff.begin(), scnt.begin(), sdsp.begin(), SData__.begin(), rcnt.begin(), rdsp.begin());
      #pragma omp parallel for
      for (sctl::Long i = 0; i < N; i++) SData__[rdsp[rank] + i] = Data[i];
    }
  }

  TRglb.ReInit(0);
  SRglb.ReInit(0);
  TSPair.ReInit(0);

    // use TData_ and SData__ only in this stage
    // step 1 build hash multimap, using <mortonid,index> as <key,value> pair
    const int threadSize=omp_get_max_threads();

    // step 2 for each target loop for its neighbor
  // fprintf(stderr,"pair-building loop in nearinterac %d\n",TData_.Dim());
    const auto TSize=TData_.Dim();
    #pragma omp parallel
    {
      const int threadId = omp_get_thread_num();
      const int threadNumber = omp_get_num_threads();
      auto TSPairThread = TSPair; // fill each local TSPair for each thread
      //WARNING: in this loop several threads may work on the same target index
      // the range is [lb,ub)
      const auto lb = (threadId + 0) * TSize / threadNumber;
      const auto ub = std::min((threadId + 1) * TSize / threadNumber,TSize);

      for (sctl::Long t0 = lb; t0 < ub;) {
        Long t1 = t0;
        MID Tmid = TData_[t0].mid;
        sctl::Vector<MID> nbrlst;
        NbrList(nbrlst, Tmid, depth, (sctl::ConstIterator<Real>)period_length0);
        while (t1 < ub && TData_[t1].mid == Tmid) t1++;
      
        for (sctl::Integer j = 0; j < nbrlst.Dim(); j++) {
          ObjData src_key;
          src_key.mid = nbrlst[j];

          // search in source
          Long s0 = std::lower_bound(SData__.begin(), SData__.end(), src_key) - SData__.begin();
          Long s1 = s0;
          if (s0 < SData__.Dim() && SData__[s0].mid == src_key.mid) {
            while (s1 < SData__.Dim() && SData__[s1].mid == src_key.mid) s1++;
          }

            for (sctl::Long t = t0; t < t1; t++) {
              for (sctl::Long s = s0; s < s1; s++) {
                Real r2 = 0, Rnear = SData__[s].rad + TData_[t].rad;
                for (sctl::Integer k = 0; k < DIM; k++) {
                  Real dx = SData__[s].coord[k] - TData_[t].coord[k];
                  if (period_length[k] > 0 && fabs(dx) > 0.5 * period_length[k]) {
                    dx -= round(dx / period_length[k]) * period_length[k];
                  }
                  r2 += dx * dx;
                }
                if (r2 > 0 && r2 < Rnear * Rnear) {
                  Long trg_idx = TData_[t].Rglb;
                  Long src_idx = SData__[s].Rglb;
                  TSPairThread.PushBack(std::pair<Long,Long>(trg_idx,src_idx));
                }
              }
            }
          }
          t0 = t1;
        }
        #pragma omp critical
        { // copy to the unified TSPair
          for(auto &pair:TSPairThread){
            TSPair.PushBack(pair);
          }
        }
      }

  // fprintf(stderr,"pair-building loop complete in nearinterac\n");
  if(TSPair.Dim()==0){
    // no pair detected
    return;
  }

  // builds SRglb and TRglb
  // fprintf(stderr,"pair sorting in nearinterac\n");
  sctl::omp_par::merge_sort(TSPair.begin(), TSPair.end());
  // fprintf(stderr,"pair sorting complete in nearinterac\n");

  #pragma omp sections
  {
    #pragma omp section
    {// one thread for TRglb
      auto lastValue=TSPair[0].first;
      TRglb.PushBack(lastValue);
      for(sctl::Long i=0;i<TSPair.Dim();i++){
        if(TSPair[i].first!=lastValue){
          lastValue=TSPair[i].first;
          TRglb.PushBack(lastValue);
        }
      }
    }

    #pragma omp section
    {// one thread for SRglb
      auto lastValue=TSPair[0].second;
      SRglb.PushBack(lastValue);
      for(sctl::Long i=0;i<TSPair.Dim();i++){
        if(TSPair[i].second!=lastValue){
          lastValue=TSPair[i].second;
          SRglb.PushBack(lastValue);
        }
      }
    }
  }
  sctl::omp_par::merge_sort(TRglb.begin(), TRglb.end());
  sctl::omp_par::merge_sort(SRglb.begin(), SRglb.end());
  // fprintf(stderr,"Rglb sorting complete in nearinterac\n");

  trg_src_pair.resize(TSPair.Dim());
  const auto TSPairSize = TSPair.Dim();
  #pragma omp parallel for
  for (sctl::Long i = 0; i < TSPairSize; i++) { // Set trg_src_pair
    const auto& TS = TSPair[i];
    auto trg_idx = TS.first;
    auto src_idx = TS.second;
    Long t = std::lower_bound(TRglb.begin(), TRglb.end(), trg_idx) - TRglb.begin();
    Long s = std::lower_bound(SRglb.begin(), SRglb.end(), src_idx) - SRglb.begin();
    SCTL_ASSERT(t < TRglb.Dim() && TRglb[t] == trg_idx);
    SCTL_ASSERT(s < SRglb.Dim() && SRglb[s] == src_idx);
    trg_src_pair[i] = std::pair<Long, Long>(t, s);

    //Real r2 = 0;
    //for (sctl::Integer k = 0; k < DIM; k++) r2 += (Svec[s].Coord()[k] - Tvec[t].Coord()[k]) * (Svec[s].Coord()[k] - Tvec[t].Coord()[k]);
    //SCTL_ASSERT(sqrt(r2) <= Svec[s].Rad() + Tvec[t].Rad());
  }

  // fprintf(stderr,"setup near complete\n");

}

template <class Real, int DIM> template <class ObjType> void NearInteraction<Real, DIM>::ForwardScatterSrc(const std::vector<ObjType>& in, std::vector<ObjType>& out) const {
  const auto& trg_src_interac = GetInteractionList();
  ForwardScatter<ObjType>(in, out, SRglb);
}

template <class Real, int DIM> template <class ObjType> void NearInteraction<Real, DIM>::ForwardScatterTrg(const std::vector<ObjType>& in, std::vector<ObjType>& out) const {
  const auto& trg_src_interac = GetInteractionList();
  ForwardScatter<ObjType>(in, out, TRglb);
}

template <class Real, int DIM> template <class ObjType> void NearInteraction<Real, DIM>::ReverseScatterTrg(const std::vector<ObjType>& in, std::vector<ObjType>& out) const {
  const auto& trg_src_interac = GetInteractionList();
  ReverseScatter<ObjType>(in, out, TRglb);
}

template <class Real, int DIM> template <class ObjType> void NearInteraction<Real, DIM>::ForwardScatter(const std::vector<ObjType>& in_vec, std::vector<ObjType>& out_vec, const sctl::Vector<Long>& recv_idx) const {
  sctl::Integer np = comm_.Size();
  sctl::Integer rank = comm_.Rank();

  sctl::Vector<Long> mins;
  { // Set mins
    mins.ReInit(np);
    Long Rglb, Nloc = in_vec.size();
    comm_.Scan<Long>(sctl::Ptr2ConstItr<Long>(&Nloc,1), sctl::Ptr2Itr<Long>(&Rglb,1), 1, sctl::Comm::CommOp::SUM);
    Rglb -= Nloc;
    comm_.Allgather<Long>(sctl::Ptr2ConstItr<Long>(&Rglb,1), 1, mins.begin(), 1);
  }

  sctl::Vector<Long> rcnt0(np), rdsp0(np);
  sctl::Vector<Long> scnt0(np), sdsp0(np);
  { // Set rcnt0, rdsp0, scnt0, sdsp0
    for (sctl::Integer i = 0; i < np; i++) {
      Long idx0 = (i + 0 < np ? std::lower_bound(recv_idx.begin(), recv_idx.end(), mins[i+0]) - recv_idx.begin() : recv_idx.Dim());
      Long idx1 = (i + 1 < np ? std::lower_bound(recv_idx.begin(), recv_idx.end(), mins[i+1]) - recv_idx.begin() : recv_idx.Dim());
      rcnt0[i] = idx1 - idx0;
      rdsp0[i] = idx0;
    }
    comm_.Alltoall(rcnt0.begin(), 1, scnt0.begin(), 1);
    sdsp0[0] = 0;
    sctl::omp_par::scan(scnt0.begin(), sdsp0.begin(), np);
  }

  sctl::Vector<Long> send_idx(sdsp0[np - 1] + scnt0[np - 1]);
  comm_.Alltoallv<Long>(recv_idx.begin(), rcnt0.begin(), rdsp0.begin(), send_idx.begin(), scnt0.begin(), sdsp0.begin());
  for (auto& idx : send_idx) { // Set send_idx <-- send_idx - mins[rank]
    idx = idx - mins[rank];
    SCTL_ASSERT(0 <= idx);
    SCTL_ASSERT(idx < in_vec.size());
  }

  sctl::Vector<char> send_buff;
  sctl::Vector<Long> ssize(send_idx.Dim());
  sctl::Vector<Long> rsize(recv_idx.Dim());
  { // Init send_buff, ssize, rsize
    std::vector<char> tmp_buff;
    for (sctl::Long i = 0; i < send_idx.Dim(); i++) { // Set send_buff, ssize
      tmp_buff.clear();
      in_vec[send_idx[i]].Pack(tmp_buff);
      for (sctl::Long j = 0; j < tmp_buff.size(); j++) send_buff.PushBack(tmp_buff[j]);
      ssize[i] = tmp_buff.size();
    }
    { // Set rsize
      comm_.Alltoallv<Long>(ssize.begin(), scnt0.begin(), sdsp0.begin(), rsize.begin(), rcnt0.begin(), rdsp0.begin());
    }
  }

  sctl::Vector<Long> scnt1(np), sdsp1(np);
  sctl::Vector<Long> rcnt1(np), rdsp1(np);
  { // Set scnt1, sdsp1
    sdsp1[0] = 0;
    for (sctl::Long i = 0; i < np; i++) {
      scnt1[i] = 0;
      for (sctl::Long j = 0; j < scnt0[i]; j++) {
        scnt1[i] += ssize[sdsp0[i] + j];
      }
    }
    sctl::omp_par::scan(scnt1.begin(), sdsp1.begin(), np);
  }
  { // Set rcnt1, rdsp1
    rdsp1[0] = 0;
    comm_.Alltoall(scnt1.begin(), 1, rcnt1.begin(), 1);
    sctl::omp_par::scan(rcnt1.begin(), rdsp1.begin(), np);
  }

  sctl::Vector<char> recv_buff(rdsp1[np - 1] + rcnt1[np - 1]);
  comm_.Alltoallv(send_buff.begin(), scnt1.begin(), sdsp1.begin(), recv_buff.begin(), rcnt1.begin(), rdsp1.begin());

  { // out_vec <-- Unpack(recv_buff)
    Long Nelem = recv_idx.Dim();

    sctl::Vector<Long> rdisp(Nelem);
    SCTL_ASSERT(rsize.Dim() == Nelem);
    if (Nelem) { // Set rdisp <-- scan(rsize)
      rdisp[0] = 0;
      sctl::omp_par::scan(rsize.begin(), rdisp.begin(), recv_idx.Dim());
    }

    out_vec.resize(Nelem);
    std::vector<char> tmp_buff;
    for (sctl::Long i = 0; i < Nelem; i++) { // Unpack recv_buff
      Long Nsize = rsize[i];
      tmp_buff.resize(Nsize);
      for (sctl::Long j = 0; j < Nsize; j++) tmp_buff[j] = recv_buff[rdisp[i] + j];
      out_vec[i].Unpack(tmp_buff);
    }
  }
}

template <class Real, int DIM> template <class ObjType> void NearInteraction<Real, DIM>::ReverseScatter(const std::vector<ObjType>& in_vec, std::vector<ObjType>& out_vec, const sctl::Vector<Long>& send_idx) const {
  SCTL_ASSERT(in_vec.size() == send_idx.Dim());
  sctl::Integer np = comm_.Size();
  sctl::Integer rank = comm_.Rank();

  sctl::Vector<Long> mins;
  { // Set mins
    mins.ReInit(np);
    Long Rglb, Nloc = out_vec.size();
    comm_.Scan<Long>(sctl::Ptr2ConstItr<Long>(&Nloc,1), sctl::Ptr2Itr<Long>(&Rglb,1), 1, sctl::Comm::CommOp::SUM);
    Rglb -= Nloc;
    comm_.Allgather<Long>(sctl::Ptr2ConstItr<Long>(&Rglb,1), 1, mins.begin(), 1);
  }

  sctl::Vector<Long> scnt0(np), sdsp0(np);
  sctl::Vector<Long> rcnt0(np), rdsp0(np);
  { // Set scnt0, sdsp0, rcnt0, rdsp0
    for (sctl::Integer i = 0; i < np; i++) {
      Long idx0 = (i + 0 < np ? std::lower_bound(send_idx.begin(), send_idx.end(), mins[i+0]) - send_idx.begin() : send_idx.Dim());
      Long idx1 = (i + 1 < np ? std::lower_bound(send_idx.begin(), send_idx.end(), mins[i+1]) - send_idx.begin() : send_idx.Dim());
      scnt0[i] = idx1 - idx0;
      sdsp0[i] = idx0;
    }
    comm_.Alltoall(scnt0.begin(), 1, rcnt0.begin(), 1);
    rdsp0[0] = 0;
    sctl::omp_par::scan(rcnt0.begin(), rdsp0.begin(), np);
  }

  sctl::Vector<Long> recv_idx(rdsp0[np - 1] + rcnt0[np - 1]);
  comm_.Alltoallv<Long>(send_idx.begin(), scnt0.begin(), sdsp0.begin(), recv_idx.begin(), rcnt0.begin(), rdsp0.begin());
  for (auto& idx : recv_idx) { // Set recv_idx <-- recv_idx - mins[rank]
    idx = idx - mins[rank];
    SCTL_ASSERT(0 <= idx);
    SCTL_ASSERT(idx < out_vec.size());
  }

  sctl::Vector<char> send_buff;
  sctl::Vector<Long> ssize(sdsp0[np - 1] + scnt0[np - 1]);
  sctl::Vector<Long> rsize(rdsp0[np - 1] + rcnt0[np - 1]);
  { // Init send_buff, ssize, rsize
    std::vector<char> tmp_buff;
    for (sctl::Long i = 0; i < in_vec.size(); i++) { // Set send_buff, ssize
      tmp_buff.clear();
      in_vec[i].Pack(tmp_buff);
      for (sctl::Long j = 0; j < tmp_buff.size(); j++) send_buff.PushBack(tmp_buff[j]);
      ssize[i] = tmp_buff.size();
    }
    { // Set rsize
      comm_.Alltoallv<Long>(ssize.begin(), scnt0.begin(), sdsp0.begin(), rsize.begin(), rcnt0.begin(), rdsp0.begin());
    }
  }

  sctl::Vector<Long> scnt1(np), sdsp1(np);
  sctl::Vector<Long> rcnt1(np), rdsp1(np);
  { // Set scnt1, sdsp1
    sdsp1[0] = 0;
    for (sctl::Long i = 0; i < np; i++) {
      scnt1[i] = 0;
      for (sctl::Long j = 0; j < scnt0[i]; j++) {
        scnt1[i] += ssize[sdsp0[i] + j];
      }
    }
    sctl::omp_par::scan(scnt1.begin(), sdsp1.begin(), np);
  }
  { // Set rcnt1, rdsp1
    rdsp1[0] = 0;
    comm_.Alltoall(scnt1.begin(), 1, rcnt1.begin(), 1);
    sctl::omp_par::scan(rcnt1.begin(), rdsp1.begin(), np);
  }

  sctl::Vector<char> recv_buff(rdsp1[np - 1] + rcnt1[np - 1]);
  comm_.Alltoallv(send_buff.begin(), scnt1.begin(), sdsp1.begin(), recv_buff.begin(), rcnt1.begin(), rdsp1.begin());

  { // out_vec <-- Unpack(recv_buff)
    Long Nelem = recv_idx.Dim();

    sctl::Vector<Long> rdisp(Nelem);
    SCTL_ASSERT(rsize.Dim() == Nelem);
    if (Nelem) { // Set rdisp <-- scan(rsize)
      rdisp[0] = 0;
      sctl::omp_par::scan(rsize.begin(), rdisp.begin(), recv_idx.Dim());
    }

    std::vector<char> tmp_buff;
    for (sctl::Long i = 0; i < Nelem; i++) { // Unpack recv_buff
      Long Nsize = rsize[i];
      tmp_buff.resize(Nsize);
      for (sctl::Long j = 0; j < Nsize; j++) tmp_buff[j] = recv_buff[rdisp[i] + j];
      out_vec[recv_idx[i]].Unpack(tmp_buff);
    }
  }
}

