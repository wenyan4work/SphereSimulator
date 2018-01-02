#ifndef INCLUDE_CUSTOMKERNEL_H_
#define INCLUDE_CUSTOMKERNEL_H_

#include <cmath>
#include <cstdlib>
#include <vector>

// pvfmm headers
#include "cheb_utils.hpp"
#include "intrin_wrapper.hpp"
#include "matrix.hpp"
#include "mem_mgr.hpp"
#include "precomp_mat.hpp"
#include "profile.hpp"
#include "vector.hpp"

namespace pvfmm {

// Stokes double layer, source: 9, target: 3
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t> >
void stokes_double_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                           Matrix<Real_t> &trg_value) {
#define SRC_BLK 500
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
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
    generic_kernel<Real_t, 9, 3, stokes_double_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t> > >(           \
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

// Stokes double layer pressure, source: 9, target: 1
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t> >
void stokes_double_pressure_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                                    Matrix<Real_t> &trg_value) {
#define SRC_BLK 500
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
        NWTN_ITER = 3;

    Real_t nwtn_scal = 1; // scaling factor for newton iterations
    for (int i = 0; i < NWTN_ITER; i++) {
        nwtn_scal = 2 * nwtn_scal * nwtn_scal * nwtn_scal;
    }
    const Real_t OOEP = 1.0 / (2 * const_pi<Real_t>());
    Vec_t inv_nwtn_scal5 = set_intrin<Vec_t, Real_t>(1.0 / (nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal));
    Vec_t negthree = set_intrin<Vec_t, Real_t>(-3.0);

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

            Vec_t tv = zero_intrin<Vec_t>();
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

                Vec_t pressure = mul_intrin(sv0, mul_intrin(dx, dx));
                pressure = add_intrin(pressure, mul_intrin(sv1, mul_intrin(dx, dy)));
                pressure = add_intrin(pressure, mul_intrin(sv2, mul_intrin(dx, dz)));
                pressure = add_intrin(pressure, mul_intrin(sv3, mul_intrin(dy, dx)));
                pressure = add_intrin(pressure, mul_intrin(sv4, mul_intrin(dy, dy)));
                pressure = add_intrin(pressure, mul_intrin(sv5, mul_intrin(dy, dz)));
                pressure = add_intrin(pressure, mul_intrin(sv6, mul_intrin(dz, dx)));
                pressure = add_intrin(pressure, mul_intrin(sv7, mul_intrin(dz, dy)));
                pressure = add_intrin(pressure, mul_intrin(sv8, mul_intrin(dz, dz)));
                pressure = mul_intrin(pressure, negthree);
                pressure = add_intrin(pressure, mul_intrin(sv0, r2));
                pressure = add_intrin(pressure, mul_intrin(sv4, r2));
                pressure = add_intrin(pressure, mul_intrin(sv8, r2));

                tv = add_intrin(tv, mul_intrin(rinv5, pressure));
            }
            Vec_t ooep = set_intrin<Vec_t, Real_t>(OOEP);

            tv = add_intrin(mul_intrin(tv, ooep), load_intrin<Vec_t>(&trg_value[0][t]));

            store_intrin(&trg_value[0][t], tv);
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
void stokes_double_pressure(T *r_src, int src_cnt, T *v_src, int dof, T *r_trg, int trg_cnt, T *v_trg,
                            mem::MemoryManager *mem_mgr) {
#define STK_KER_NWTN(nwtn)                                                                                             \
    if (newton_iter == nwtn)                                                                                           \
    generic_kernel<Real_t, 9, 1, stokes_double_pressure_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t> > >(  \
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

// Stokes double layer gradient of velocity, source: 9, target: 9
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t> >
void stokes_dgrad_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                          Matrix<Real_t> &trg_value) {
#define SRC_BLK 100
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
        NWTN_ITER = 3;

    Real_t nwtn_scal = 1; // scaling factor for newton iterations
    for (int i = 0; i < NWTN_ITER; i++) {
        nwtn_scal = 2 * nwtn_scal * nwtn_scal * nwtn_scal;
    }
    const Real_t OOEP = -3.0 / (4 * const_pi<Real_t>());
    Vec_t inv_nwtn_scal7 = set_intrin<Vec_t, Real_t>(
        1.0 / (nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal));
    Vec_t inv_nwtn_scal5 = set_intrin<Vec_t, Real_t>(1.0 / (nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal * nwtn_scal));
    Vec_t nwtn_scal2 = set_intrin<Vec_t, Real_t>((nwtn_scal * nwtn_scal));
    Vec_t negone = set_intrin<Vec_t, Real_t>(-1.0);
    Vec_t negtwo = set_intrin<Vec_t, Real_t>(-2.0);
    Vec_t five = set_intrin<Vec_t, Real_t>(5.0);

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

            Vec_t tv0ac = zero_intrin<Vec_t>();
            Vec_t tv1ac = zero_intrin<Vec_t>();
            Vec_t tv2ac = zero_intrin<Vec_t>();
            Vec_t tv3ac = zero_intrin<Vec_t>();
            Vec_t tv4ac = zero_intrin<Vec_t>();
            Vec_t tv5ac = zero_intrin<Vec_t>();
            Vec_t tv6ac = zero_intrin<Vec_t>();
            Vec_t tv7ac = zero_intrin<Vec_t>();
            Vec_t tv8ac = zero_intrin<Vec_t>();

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

                Vec_t rinv5 = mul_intrin(rinv, rinv4);
                Vec_t rinv7 = mul_intrin(mul_intrin(rinv2, rinv5), inv_nwtn_scal7);
                rinv5 = mul_intrin(rinv5, inv_nwtn_scal5);

                Vec_t commonCoeff = mul_intrin(sv0, mul_intrin(dx, dx));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv1, mul_intrin(dx, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv2, mul_intrin(dx, dz)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv3, mul_intrin(dy, dx)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv4, mul_intrin(dy, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv5, mul_intrin(dy, dz)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv6, mul_intrin(dz, dx)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv7, mul_intrin(dz, dy)));
                commonCoeff = add_intrin(commonCoeff, mul_intrin(sv8, mul_intrin(dz, dz)));
                commonCoeff = mul_intrin(negone, commonCoeff);

                // (-2*(t0 - s0)*sv0 - (t1 - s1)*(sv1 + sv3) - (t2 - s2)*(sv2 + sv6))
                Vec_t dcFd0 = mul_intrin(negtwo, mul_intrin(dx, sv0));
                dcFd0 = add_intrin(dcFd0, mul_intrin(negone, mul_intrin(dy, add_intrin(sv1, sv3))));
                dcFd0 = add_intrin(dcFd0, mul_intrin(negone, mul_intrin(dz, add_intrin(sv2, sv6))));

                // (-2*(t1 - s1)*sv4 - (t0 - s0)*(sv1 + sv3) - (t2 - s2)*(sv5 + sv7))
                Vec_t dcFd1 = mul_intrin(negtwo, mul_intrin(dy, sv4));
                dcFd1 = add_intrin(dcFd1, mul_intrin(negone, mul_intrin(dx, add_intrin(sv1, sv3))));
                dcFd1 = add_intrin(dcFd1, mul_intrin(negone, mul_intrin(dz, add_intrin(sv5, sv7))));

                // (-2*(t2 - s2)*sv8 - (t0 - s0)*(sv2 + sv6) - (t1 - s1)*(sv5 + sv7))
                Vec_t dcFd2 = mul_intrin(negtwo, mul_intrin(dz, sv8));
                dcFd2 = add_intrin(dcFd2, mul_intrin(negone, mul_intrin(dx, add_intrin(sv2, sv6))));
                dcFd2 = add_intrin(dcFd2, mul_intrin(negone, mul_intrin(dy, add_intrin(sv5, sv7))));

                Vec_t tv0 = zero_intrin<Vec_t>();
                Vec_t tv1 = zero_intrin<Vec_t>();
                Vec_t tv2 = zero_intrin<Vec_t>();
                Vec_t tv3 = zero_intrin<Vec_t>();
                Vec_t tv4 = zero_intrin<Vec_t>();
                Vec_t tv5 = zero_intrin<Vec_t>();
                Vec_t tv6 = zero_intrin<Vec_t>();
                Vec_t tv7 = zero_intrin<Vec_t>();
                Vec_t tv8 = zero_intrin<Vec_t>();

                // (5 * rrtensor * commonCoeff / rnorm ^ 7
                //  - Outer[Times, rvec, {dcFd0, dcFd1, dcFd2}] / rnorm ^5
                //  - IdentityMatrix[3] * commonCoeff / rnorm ^ 5);
                tv0 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dx, dx)));
                tv0 = mul_intrin(tv0, rinv7);
                tv0 = add_intrin(tv0, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dx, dcFd0))));
                tv0 = add_intrin(tv0, mul_intrin(negone, mul_intrin(rinv5, commonCoeff)));
                tv1 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dx, dy)));
                tv1 = mul_intrin(tv1, rinv7);
                tv1 = add_intrin(tv1, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dx, dcFd1))));
                tv2 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dx, dz)));
                tv2 = mul_intrin(tv2, rinv7);
                tv2 = add_intrin(tv2, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dx, dcFd2))));

                tv3 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dy, dx)));
                tv3 = mul_intrin(tv3, rinv7);
                tv3 = add_intrin(tv3, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dy, dcFd0))));
                tv4 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dy, dy)));
                tv4 = mul_intrin(tv4, rinv7);
                tv4 = add_intrin(tv4, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dy, dcFd1))));
                tv4 = add_intrin(tv4, mul_intrin(negone, mul_intrin(rinv5, commonCoeff)));
                tv5 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dy, dz)));
                tv5 = mul_intrin(tv5, rinv7);
                tv5 = add_intrin(tv5, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dy, dcFd2))));

                tv6 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dz, dx)));
                tv6 = mul_intrin(tv6, rinv7);
                tv6 = add_intrin(tv6, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dz, dcFd0))));
                tv7 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dz, dy)));
                tv7 = mul_intrin(tv7, rinv7);
                tv7 = add_intrin(tv7, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dz, dcFd1))));
                tv8 = mul_intrin(five, mul_intrin(commonCoeff, mul_intrin(dz, dz)));
                tv8 = mul_intrin(tv8, rinv7);
                tv8 = add_intrin(tv8, mul_intrin(negone, mul_intrin(rinv5, mul_intrin(dz, dcFd2))));
                tv8 = add_intrin(tv8, mul_intrin(negone, mul_intrin(rinv5, commonCoeff)));

                tv0ac = add_intrin(tv0ac, tv0);
                tv1ac = add_intrin(tv1ac, tv1);
                tv2ac = add_intrin(tv2ac, tv2);
                tv3ac = add_intrin(tv3ac, tv3);
                tv4ac = add_intrin(tv4ac, tv4);
                tv5ac = add_intrin(tv5ac, tv5);
                tv6ac = add_intrin(tv6ac, tv6);
                tv7ac = add_intrin(tv7ac, tv7);
                tv8ac = add_intrin(tv8ac, tv8);
            }
            Vec_t ooep = set_intrin<Vec_t, Real_t>(OOEP);

            tv0ac = add_intrin(mul_intrin(tv0ac, ooep), load_intrin<Vec_t>(&trg_value[0][t]));
            tv1ac = add_intrin(mul_intrin(tv1ac, ooep), load_intrin<Vec_t>(&trg_value[1][t]));
            tv2ac = add_intrin(mul_intrin(tv2ac, ooep), load_intrin<Vec_t>(&trg_value[2][t]));
            tv3ac = add_intrin(mul_intrin(tv3ac, ooep), load_intrin<Vec_t>(&trg_value[3][t]));
            tv4ac = add_intrin(mul_intrin(tv4ac, ooep), load_intrin<Vec_t>(&trg_value[4][t]));
            tv5ac = add_intrin(mul_intrin(tv5ac, ooep), load_intrin<Vec_t>(&trg_value[5][t]));
            tv6ac = add_intrin(mul_intrin(tv6ac, ooep), load_intrin<Vec_t>(&trg_value[6][t]));
            tv7ac = add_intrin(mul_intrin(tv7ac, ooep), load_intrin<Vec_t>(&trg_value[7][t]));
            tv8ac = add_intrin(mul_intrin(tv8ac, ooep), load_intrin<Vec_t>(&trg_value[8][t]));

            store_intrin(&trg_value[0][t], tv0ac);
            store_intrin(&trg_value[1][t], tv1ac);
            store_intrin(&trg_value[2][t], tv2ac);
            store_intrin(&trg_value[3][t], tv3ac);
            store_intrin(&trg_value[4][t], tv4ac);
            store_intrin(&trg_value[5][t], tv5ac);
            store_intrin(&trg_value[6][t], tv6ac);
            store_intrin(&trg_value[7][t], tv7ac);
            store_intrin(&trg_value[8][t], tv8ac);
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
void stokes_dgrad(T *r_src, int src_cnt, T *v_src, int dof, T *r_trg, int trg_cnt, T *v_trg,
                  mem::MemoryManager *mem_mgr) {
#define STK_KER_NWTN(nwtn)                                                                                             \
    if (newton_iter == nwtn)                                                                                           \
    generic_kernel<Real_t, 9, 9, stokes_dgrad_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t> > >(            \
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
template <class Real_t, class Vec_t = Real_t, Vec_t (*RSQRT_INTRIN)(Vec_t) = rsqrt_intrin0<Vec_t> >
void stokes_traction_uKernel(Matrix<Real_t> &src_coord, Matrix<Real_t> &src_value, Matrix<Real_t> &trg_coord,
                             Matrix<Real_t> &trg_value) {
#define SRC_BLK 500
    size_t VecLen = sizeof(Vec_t) / sizeof(Real_t);

    //// Number of newton iterations
    size_t NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin0<Vec_t, Real_t>)
        NWTN_ITER = 0;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin1<Vec_t, Real_t>)
        NWTN_ITER = 1;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin2<Vec_t, Real_t>)
        NWTN_ITER = 2;
    if (RSQRT_INTRIN == (Vec_t (*)(Vec_t))rsqrt_intrin3<Vec_t, Real_t>)
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

                tv0 = add_intrin(tv0, mul_intrin(rinv5, mul_intrin(mul_intrin(dx, dx), commonCoeff)));
                tv1 = add_intrin(tv1, mul_intrin(rinv5, mul_intrin(mul_intrin(dx, dy), commonCoeff)));
                tv2 = add_intrin(tv2, mul_intrin(rinv5, mul_intrin(mul_intrin(dx, dz), commonCoeff)));
                tv3 = add_intrin(tv3, mul_intrin(rinv5, mul_intrin(mul_intrin(dy, dx), commonCoeff)));
                tv4 = add_intrin(tv4, mul_intrin(rinv5, mul_intrin(mul_intrin(dy, dy), commonCoeff)));
                tv5 = add_intrin(tv5, mul_intrin(rinv5, mul_intrin(mul_intrin(dy, dz), commonCoeff)));
                tv6 = add_intrin(tv6, mul_intrin(rinv5, mul_intrin(mul_intrin(dz, dx), commonCoeff)));
                tv7 = add_intrin(tv7, mul_intrin(rinv5, mul_intrin(mul_intrin(dz, dy), commonCoeff)));
                tv8 = add_intrin(tv8, mul_intrin(rinv5, mul_intrin(mul_intrin(dz, dz), commonCoeff)));
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
    generic_kernel<Real_t, 3, 9, stokes_traction_uKernel<Real_t, Vec_t, rsqrt_intrin##nwtn<Vec_t, Real_t> > >(         \
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
    inline static const Kernel<T> &Traction();
    inline static const Kernel<T> &Double();
    inline static const Kernel<T> &DoublePressure();
    inline static const Kernel<T> &DoubleGrad();
};

// reference, how to define a kernel only for S2T/M2T/L2T
// template <class T>
// const Kernel<T> &LaplaceCustomKernel<T>::potentialGradient() {
//    static Kernel<T> potn_ker = BuildKernel<T, laplace_poten<T, 1> >("laplace", 3, std::pair<int, int>(1, 1));
//    static Kernel<T> pgrad_ker =
//        BuildKernel<T, laplace_grad<T, 1> >("laplace_pgrad", 3, std::pair<int, int>(1, 4), &potn_ker, &potn_ker, NULL,
//                                            &potn_ker, &potn_ker, NULL, &potn_ker, NULL);
//    return pgrad_ker;
//}

#define NEWTON_ITE 1
template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::Traction() {
    static Kernel<T> stk_ker = BuildKernel<T, stokes_vel<T, NEWTON_ITE> >("stokes_vel", 3, std::pair<int, int>(3, 3));
    static Kernel<T> trac_ker =
        BuildKernel<T, stokes_traction<T, NEWTON_ITE> >("stokes_traction", 3, std::pair<int, int>(3, 9), &stk_ker,
                                                        &stk_ker, NULL, &stk_ker, &stk_ker, NULL, &stk_ker, NULL);
    return trac_ker;
}

template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::Double() {
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    return double_ker;
}

template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::DoublePressure() {
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    static Kernel<T> pressure_ker = BuildKernel<T, stokes_double_pressure<T, NEWTON_ITE> >(
        "stokes_dpressure", 3, std::pair<int, int>(9, 1), &double_ker, &double_ker, NULL, &double_ker, &double_ker,
        NULL, &double_ker, NULL);
    return pressure_ker;
}

template <class T>
inline const Kernel<T> &StokesCustomKernel<T>::DoubleGrad() {
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    static Kernel<T> grad_ker = BuildKernel<T, stokes_dgrad<T, NEWTON_ITE> >(
        "stokes_dgrad", 3, std::pair<int, int>(9, 9), &double_ker, &double_ker, NULL, &double_ker, &double_ker, NULL,
        &double_ker, NULL);
    return grad_ker;
}
#undef NEWTON_ITE

// Two Newton Iterations for double precision
#define NEWTON_ITE 2
template <>
inline const Kernel<double> &StokesCustomKernel<double>::Traction() {
    using T = double;
    static Kernel<T> stk_ker = BuildKernel<T, stokes_vel<T, NEWTON_ITE> >("stokes_vel", 3, std::pair<int, int>(3, 3));
    static Kernel<T> trac_ker =
        BuildKernel<T, stokes_traction<T, NEWTON_ITE> >("stokes_traction", 3, std::pair<int, int>(3, 9), &stk_ker,
                                                        &stk_ker, NULL, &stk_ker, &stk_ker, NULL, &stk_ker, NULL);
    return trac_ker;
}

template <>
inline const Kernel<double> &StokesCustomKernel<double>::Double() {
    using T = double;
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    return double_ker;
}

template <>
inline const Kernel<double> &StokesCustomKernel<double>::DoublePressure() {
    using T = double;
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    static Kernel<T> pressure_ker = BuildKernel<T, stokes_double_pressure<T, NEWTON_ITE> >(
        "stokes_dpressure", 3, std::pair<int, int>(9, 1), &double_ker, &double_ker, NULL, &double_ker, &double_ker,
        NULL, &double_ker, NULL);
    return pressure_ker;
}

template <>
inline const Kernel<double> &StokesCustomKernel<double>::DoubleGrad() {
    using T = double;
    static Kernel<T> double_ker =
        BuildKernel<T, stokes_double<T, NEWTON_ITE> >("stokes_double", 3, std::pair<int, int>(9, 3));
    static Kernel<T> grad_ker = BuildKernel<T, stokes_dgrad<T, NEWTON_ITE> >(
        "stokes_dgrad", 3, std::pair<int, int>(9, 9), &double_ker, &double_ker, NULL, &double_ker, &double_ker, NULL,
        &double_ker, NULL);
    return grad_ker;
}
#undef NEWTON_ITE
}
#endif
