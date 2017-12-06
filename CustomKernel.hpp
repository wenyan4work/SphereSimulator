#ifndef INCLUDE_CUSTOMKERNEL_H_
#define INCLUDE_CUSTOMKERNEL_H_

#include <cmath>
#include <cstdlib>
#include <vector>

// pvfmm headers
#include <cheb_utils.hpp>
#include <intrin_wrapper.hpp>
#include <matrix.hpp>
#include <mem_mgr.hpp>
#include <precomp_mat.hpp>
#include <profile.hpp>
#include <vector.hpp>

namespace pvfmm {

// Stokes double layer, source: 9, target: 3
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t>>
void stokes_double_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                           Matrix<Real_t> &trg_value) {
#define SRC_BLK 500
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
        NWTN_ITER = 3;

    Real_t nwtn_scal = 1; // scaling factor for newton iterations
    for (int i = 0; i < NWTN_ITER; i++) {
        nwtn_scal = 2 * nwtn_scal * nwtn_scal * nwtn_scal;
    }
    const Real_t OOEP = -3.0 / (4 * const_pi<Real_t>());
    Vec_t inv_nwtn_scal5 = set_intrin<Vec_t, Real_t>(1.0 / (nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal));

    size_t src_cnt_ = src_coord.Dim(1);
    size_t trg_cnt_ = trg_coord.Dim(1);
    for (size_t sblk = 0; sblk < src_cnt_; sblk += SRC_BLK) {
        size_t src_cnt = src_cnt_ - sblk;
        if (src_cnt > SRC_BLK)
            src_cnt = SRC_BLK;
        for (size_t t = 0; t < trg_cnt_; t += VecLen) {
            Vec_t tx = load_intrin<Vec_t>(&trg_coord[0][t]);
            Vec_t ty = load_intrin<Vec_t>(&trg_coord[1][t]);
            Vec_t tz = load_intrin<Vec_t>(&trg_coord[2][t]);

            Vec_t tvx = zero_intrin<Vec_t>();
            Vec_t tvy = zero_intrin<Vec_t>();
            Vec_t tvz = zero_intrin<Vec_t>();
            for (size_t s = sblk; s < sblk + src_cnt; s++) {
                Vec_t dx = sub_intrin(tx, bcast_intrin<Vec_t>(&src_coord[0][s]));
                Vec_t dy = sub_intrin(ty, bcast_intrin<Vec_t>(&src_coord[1][s]));
                Vec_t dz = sub_intrin(tz, bcast_intrin<Vec_t>(&src_coord[2][s]));

                Vec_t sv0 = bcast_intrin<Vec_t>(&src_value[0][s]);
                Vec_t sv1 = bcast_intrin<Vec_t>(&src_value[1][s]);
                Vec_t sv2 = bcast_intrin<Vec_t>(&src_value[2][s]);
                Vec_t sv3 = bcast_intrin<Vec_t>(&src_value[3][s]);
                Vec_t sv4 = bcast_intrin<Vec_t>(&src_value[4][s]);
                Vec_t sv5 = bcast_intrin<Vec_t>(&src_value[5][s]);
                Vec_t sv6 = bcast_intrin<Vec_t>(&src_value[6][s]);
                Vec_t sv7 = bcast_intrin<Vec_t>(&src_value[7][s]);
                Vec_t sv8 = bcast_intrin<Vec_t>(&src_value[8][s]);

                Vec_t r2 = mul_intrin(dx, dx);
                r2 = add_intrin(r2, mul_intrin(dy, dy));
                r2 = add_intrin(r2, mul_intrin(dz, dz));

                Vec_t rinv = RSQRT_INTRIN(r2);
                Vec_t rinv2 = mul_intrin(rinv, rinv);
                Vec_t rinv4 = mul_intrin(rinv2, rinv2);

                Vec_t rinv5 = mul_intrin(mul_intrin(rinv, rinv4), inv_nwtn_scal5);

                Vec_t commonCoeff = mul_intrin(sv0, mul_intrin(dx, dx));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv1, mul_intrin(dx, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv2, mul_intrin(dx, dz)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv3, mul_intrin(dy, dx)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv4, mul_intrin(dy, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv5, mul_intrin(dy, dz)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv6, mul_intrin(dz, dx)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv7, mul_intrin(dz, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv8, mul_intrin(dz, dz)));

                tvx = add_intrin(tvx, mul_intrin(rinv5, mul_intrin(dx, commonCoeff)));
                tvy = add_intrin(tvy, mul_intrin(rinv5, mul_intrin(dy, commonCoeff)));
                tvz = add_intrin(tvz, mul_intrin(rinv5, mul_intrin(dz, commonCoeff)));
            }
            Vec_t ooep = set_intrin<Vec_t, Real_t>(OOEP);

            tvx = add_intrin(mul_intrin(tvx, ooep), load_intrin<Vec_t>(&trg_value[0][t]));
            tvy = add_intrin(mul_intrin(tvy, ooep), load_intrin<Vec_t>(&trg_value[1][t]));
            tvz = add_intrin(mul_intrin(tvz, ooep), load_intrin<Vec_t>(&trg_value[2][t]));

            store_intrin(&trg_value[0][t], tvx);
            store_intrin(&trg_value[1][t], tvy);
            store_intrin(&trg_value[2][t], tvz);
        }
    }

    { // Add FLOPS
#ifndef __MIC__
        Profile::Add_FLOP((long long)trg_cnt_ * (long long)src_cnt_ * (29 + 4 * (NWTN_ITER)));
#endif
    }
#undef SRC_BLK
}

template <class T, int newton_iter = 0>
void stokes_double(T *r_src, int src_cnt, T *v_src, int dof, T *r_trg, int trg_cnt, T *v_trg,
                   mem::MemoryManager *mem_mgr) {
#define STK_KER_NWTN(nwtn)                                                                                             \
    if (newton_iter == nwtn)                                                                                           \
    generic_kernel<Real_t, 9, 3, stokes_double_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t>>>(             \
        (Real_t *)r_src, src_cnt, (Real_t *)v_src, dof, (Real_t *)r_trg, trg_cnt, (Real_t *)v_trg, mem_mgr)
#define STOKES_KERNEL                                                                                                  \
    STK_KER_NWTN(0);                                                                                                   \
    STK_KER_NWTN(1);                                                                                                   \
    STK_KER_NWTN(2);                                                                                                   \
    STK_KER_NWTN(3);

    if (mem::TypeTraits<T>::ID() == mem::TypeTraits<float>::ID()) {
        typedef float Real_t;
#if defined __MIC__
#define Vec_t Real_t
#elif defined __AVX__
#define Vec_t __m256
#elif defined __SSE3__
#define Vec_t __m128
#else
#define Vec_t Real_t
#endif
        STOKES_KERNEL;
#undef Vec_t
    } else if (mem::TypeTraits<T>::ID() == mem::TypeTraits<double>::ID()) {
        typedef double Real_t;
#if defined __MIC__
#define Vec_t Real_t
#elif defined __AVX__
#define Vec_t __m256d
#elif defined __SSE3__
#define Vec_t __m128d
#else
#define Vec_t Real_t
#endif
        STOKES_KERNEL;
#undef Vec_t
    } else {
        typedef T Real_t;
#define Vec_t Real_t
        STOKES_KERNEL;
#undef Vec_t
    }

#undef STK_KER_NWTN
#undef STOKES_KERNEL
}

// Stokes traction kernel, source: 3, target: 9
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t>>
void stokes_traction_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                             Matrix<Real_t> &trg_value) {
#define SRC_BLK 500
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t(*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
        NWTN_ITER = 3;

    Real_t nwtn_scal = 1; // scaling factor for newton iterations
    for (int i = 0; i < NWTN_ITER; i++) {
        nwtn_scal = 2 * nwtn_scal * nwtn_scal * nwtn_scal;
    }
    const Real_t OOEP = -3.0 / (4 * const_pi<Real_t>());
    Vec_t inv_nwtn_scal5 = set_intrin<Vec_t, Real_t>(1.0 / (nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal));

    size_t src_cnt_ = src_coord.Dim(1);
    size_t trg_cnt_ = trg_coord.Dim(1);
    for (size_t sblk = 0; sblk < src_cnt_; sblk += SRC_BLK) {
        size_t src_cnt = src_cnt_ - sblk;
        if (src_cnt > SRC_BLK)
            src_cnt = SRC_BLK;
        for (size_t t = 0; t < trg_cnt_; t += VecLen) {
            Vec_t tx = load_intrin<Vec_t>(&trg_coord[0][t]);
            Vec_t ty = load_intrin<Vec_t>(&trg_coord[1][t]);
            Vec_t tz = load_intrin<Vec_t>(&trg_coord[2][t]);

            Vec_t tv0 = zero_intrin<Vec_t>();
            Vec_t tv1 = zero_intrin<Vec_t>();
            Vec_t tv2 = zero_intrin<Vec_t>();
            Vec_t tv3 = zero_intrin<Vec_t>();
            Vec_t tv4 = zero_intrin<Vec_t>();
            Vec_t tv5 = zero_intrin<Vec_t>();
            Vec_t tv6 = zero_intrin<Vec_t>();
            Vec_t tv7 = zero_intrin<Vec_t>();
            Vec_t tv8 = zero_intrin<Vec_t>();

            for (size_t s = sblk; s < sblk + src_cnt; s++) {
                Vec_t dx = sub_intrin(tx, bcast_intrin<Vec_t>(&src_coord[0][s]));
                Vec_t dy = sub_intrin(ty, bcast_intrin<Vec_t>(&src_coord[1][s]));
                Vec_t dz = sub_intrin(tz, bcast_intrin<Vec_t>(&src_coord[2][s]));

                Vec_t sv0 = bcast_intrin<Vec_t>(&src_value[0][s]);
                Vec_t sv1 = bcast_intrin<Vec_t>(&src_value[1][s]);
                Vec_t sv2 = bcast_intrin<Vec_t>(&src_value[2][s]);

                Vec_t r2 = mul_intrin(dx, dx);
                r2 = add_intrin(r2, mul_intrin(dy, dy));
                r2 = add_intrin(r2, mul_intrin(dz, dz));

                Vec_t rinv = RSQRT_INTRIN(r2);
                Vec_t rinv2 = mul_intrin(rinv, rinv);
                Vec_t rinv4 = mul_intrin(rinv2, rinv2);

                Vec_t rinv5 = mul_intrin(mul_intrin(rinv, rinv4), inv_nwtn_scal5);

                Vec_t commonCoeff = mul_intrin(sv0, dx);
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv1, dy));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv2, dz));

                tv0 = add_intrin(tv0, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dx, dx))));
                tv1 = add_intrin(tv1, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dx, dy))));
                tv2 = add_intrin(tv2, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dx, dz))));
                tv3 = add_intrin(tv3, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dy, dx))));
                tv4 = add_intrin(tv4, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dy, dy))));
                tv5 = add_intrin(tv5, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dy, dz))));
                tv6 = add_intrin(tv6, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dz, dx))));
                tv7 = add_intrin(tv7, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dz, dy))));
                tv8 = add_intrin(tv8, mul_intrin(rinv5, mul_intrin(commonCoeff, mul_intrin(dz, dz))));
            }
            Vec_t ooep = set_intrin<Vec_t, Real_t>(OOEP);

            tv0 = add_intrin(mul_intrin(tv0, ooep), load_intrin<Vec_t>(&trg_value[0][t]));
            tv1 = add_intrin(mul_intrin(tv1, ooep), load_intrin<Vec_t>(&trg_value[1][t]));
            tv2 = add_intrin(mul_intrin(tv2, ooep), load_intrin<Vec_t>(&trg_value[2][t]));
            tv3 = add_intrin(mul_intrin(tv3, ooep), load_intrin<Vec_t>(&trg_value[3][t]));
            tv4 = add_intrin(mul_intrin(tv4, ooep), load_intrin<Vec_t>(&trg_value[4][t]));
            tv5 = add_intrin(mul_intrin(tv5, ooep), load_intrin<Vec_t>(&trg_value[5][t]));
            tv6 = add_intrin(mul_intrin(tv6, ooep), load_intrin<Vec_t>(&trg_value[6][t]));
            tv7 = add_intrin(mul_intrin(tv7, ooep), load_intrin<Vec_t>(&trg_value[7][t]));
            tv8 = add_intrin(mul_intrin(tv8, ooep), load_intrin<Vec_t>(&trg_value[8][t]));

            store_intrin(&trg_value[0][t], tv0);
            store_intrin(&trg_value[1][t], tv1);
            store_intrin(&trg_value[2][t], tv2);
            store_intrin(&trg_value[3][t], tv3);
            store_intrin(&trg_value[4][t], tv4);
            store_intrin(&trg_value[5][t], tv5);
            store_intrin(&trg_value[6][t], tv6);
            store_intrin(&trg_value[7][t], tv7);
            store_intrin(&trg_value[8][t], tv8);
        }
    }

    { // Add FLOPS
#ifndef __MIC__
        Profile::Add_FLOP((long long)trg_cnt_ * (long long)src_cnt_ * (29 + 4 * (NWTN_ITER)));
#endif
    }
#undef SRC_BLK
}

template <class T, int newton_iter = 0>
void stokes_traction(T *r_src, int src_cnt, T *v_src, int dof, T *r_trg, int trg_cnt, T *v_trg,
                     mem::MemoryManager *mem_mgr) {
#define STK_KER_NWTN(nwtn)                                                                                             \
    if (newton_iter == nwtn)                                                                                           \
    generic_kernel<Real_t, 3, 9, stokes_traction_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t>>>(           \
        (Real_t *)r_src, src_cnt, (Real_t *)v_src, dof, (Real_t *)r_trg, trg_cnt, (Real_t *)v_trg, mem_mgr)
#define STOKES_KERNEL                                                                                                  \
    STK_KER_NWTN(0);                                                                                                   \
    STK_KER_NWTN(1);                                                                                                   \
    STK_KER_NWTN(2);                                                                                                   \
    STK_KER_NWTN(3);

    if (mem::TypeTraits<T>::ID() == mem::TypeTraits<float>::ID()) {
        typedef float Real_t;
#if defined __MIC__
#define Vec_t Real_t
#elif defined __AVX__
#define Vec_t __m256
#elif defined __SSE3__
#define Vec_t __m128
#else
#define Vec_t Real_t
#endif
        STOKES_KERNEL;
#undef Vec_t
    } else if (mem::TypeTraits<T>::ID() == mem::TypeTraits<double>::ID()) {
        typedef double Real_t;
#if defined __MIC__
#define Vec_t Real_t
#elif defined __AVX__
#define Vec_t __m256d
#elif defined __SSE3__
#define Vec_t __m128d
#else
#define Vec_t Real_t
#endif
        STOKES_KERNEL;
#undef Vec_t
    } else {
        typedef T Real_t;
#define Vec_t Real_t
        STOKES_KERNEL;
#undef Vec_t
    }

#undef STK_KER_NWTN
#undef STOKES_KERNEL
}

template <class T>
struct StokesCustomKernel {
    inline static const Kernel<T> &Double();
    inline static const Kernel<T> &Traction();
};

template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::Double() {
    static Kernel<T> potn_ker = BuildKernel<T, stokes_double<T, 1>>("stokes_double", 3, std::pair<int, int>(9, 3));
    return potn_ker;
}

template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::Traction() {
    static Kernel<T> potn_ker = BuildKernel<T, stokes_traction<T, 1>>("stokes_traction", 3, std::pair<int, int>(3, 9));
    return potn_ker;
}

// two newton iterations for double type
template <>
inline const Kernel<double> &StokesCustomKernel<double>::Double() {
    typedef double T;
    static Kernel<T> potn_ker = BuildKernel<T, stokes_double<T, 2>>("stokes_double", 3, std::pair<int, int>(9, 3));
    return potn_ker;
}

template <>
inline const Kernel<double> &StokesCustomKernel<double>::Traction() {
    typedef double T;
    static Kernel<T> potn_ker = BuildKernel<T, stokes_traction<T, 2>>("stokes_traction", 3, std::pair<int, int>(3, 9));
    return potn_ker;
}
}

#endif