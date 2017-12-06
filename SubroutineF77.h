/*
 * ForRoutines.h
 *
 *  Created on: Nov 27, 2017
 *      Author: wyan
 */

#ifndef SPHEREVWX_FORROUTINES_H_
#define SPHEREVWX_FORROUTINES_H_

// complex class
struct complx {
    double dr;
    double di;
};

/*
 * https://software.intel.com/en-us/node/678430
 * Do NOT use long double with mixing c and fortran:
 * gcc* on 32-bit Linux*: long double is an 80-bit floating type, as supported by the X87 instruction set.
 *                        Intel® Fortran does not support this, so C_LONG_DOUBLE is -1.
 * gcc on OS X*:          C_LONG_DOUBLE is defined as a 128-bit type that is the same as the Intel® Fortran data type REAL(16).
 * 64-bit Linux* and Windows*:
 *                        long double is treated the same as double, so C_LONG_DOUBLE is 8.
 *                        To ensure matching, make sure you use the constants for kind values and the corresponding types in C.
 * */

// sphtrans Fortran function declarations
extern "C" {
    //Spherical harmonic function evaluation
    void ynmeva_(int *n, double *x, double *y, double *z, struct complx Y[]);
    void ynmeva2_(int *n, double *x, double *y, double *z, struct complx Y[]);
    void xnmeva2_(int *n, double *x, double *y, double *z, struct complx X1[], struct complx X2[]);
    void unmeva2_(int *n, double *x, double *y, double *z, struct complx U1[], struct complx U2[]);
    void vnmeva2_(int *n, double *x, double *y, double *z, struct complx V1[], struct complx Y[], struct complx V2[]);
    void wnmeva2_(int *n, double *x, double *y, double *z, struct complx W1[], struct complx Y[], struct complx W2[]);

    // Legendre functions
    void legeexps_(int *i, int *n, double xs[], double u[], double v[], double ws[]);
    void legewhts_(int *n, double xs[], double ws[], bool ifwh);

    //Scalar Sph transformations
    //void mpoleinit(int *n, struct complx mpole[]);
    void sphtrans_cmpl_lege_init_(int *n, int *nphi, int *nth, double cth[], double wh[], double ynms[],
            struct complx ws[]);
    void sphtrans_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid[], double cth[],
            double ynms[], struct complx ws[]);
    void sphtrans_cmpl_lege_brute_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid[]);
    void sphtrans_fwd_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid2[], double cth[],
            double wh[], double ynms[], struct complx ws[]);
    void fftnext235_(int *n, int *n2);

    //Vector Sph transformations
    void sphtrans_xu_cmpl_lege_init_(int *n, int *nphi, int *nth, double cth[], double wh[], double ynms[],
            double dnms[], struct complx ws[]);
    void sphtrans_d_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid[], double cth[],
            double dnms[], struct complx ws[]);
    void sphtrans_x_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[], double cth[], double wh[], double ynms[], struct complx ws[]);
    void sphtrans_u_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[], double cth[], double wh[], double ynms[], struct complx ws[]);
    void sphtrans_y_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid[], double cth[],
            double ynms[], struct complx ws[]);
    void sphtrans_x_fwd_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[], double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]);
    void sphtrans_u_fwd_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[], double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]);
    void sphtrans_y_fwd_cmpl_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid[], double cth[],
            double wh[], double ynms[], struct complx ws[]);
    void sphtrans_x_cmpl_lege_brute_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[]);
    void sphtrans_u_cmpl_lege_brute_(int *n, struct complx mpole[], int *nphi, int *nth, struct complx fgrid1[],
            struct complx fgrid2[]);

    // STFMM3DLIB functions
    void stfmm3dpartself(int *ier, int *iprec, int *nsrc, double src[], bool *ifsingle, double sigma_sl[],
            bool *ifdouble, double sigma_dl[], double sigma_dv[], bool *ifpot, double pot[], double pre[], bool *ifgrad,
            double grad[]);
    void stfmm3dparttarg(int *ier, int *iprec, int *nsrc, double src[], bool *ifsingle, double sigma_sl[],
            bool *ifdouble, double sigma_dl[], double sigma_dv[], bool *ifpot, double pot[], double pre[], bool *ifgrad,
            double grad[], int *ntrg, double trg[], bool *ifpottrg, double pottrg[], double pretrg[], bool *ifgradtrg,
            double gradtrg[]);
    void st3dpartdirect(int *nsrc, double src[], bool *ifsingle, double sigma_sl[], bool *ifdouble, double sigma_dl[],
            double sigma_dv[], bool *ifpot, double pot[], double pre[], bool *ifgrad, double grad[], int *ntrg,
            double trg[], bool *ifpottrg, double pottrg[], double pretrg[], bool *ifgradtrg, double gradtrg[]);
}

#endif /* SPHEREVWX_FORROUTINES_H_ */
