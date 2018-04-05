#ifndef TRACTIONEXTERIOR_HPP_
#define TRACTIONEXTERIOR_HPP_

#include <cmath>
#include <complex>

template <class Real>
inline Real Power(const Real &a, const int &b) {
    return std::pow(a, b);
}

// this is a purely real function, can be called separately using either
// the real or the imaginary part of the spectral coefficients of phiVnm, phiWnm, phiXnm
// input: n,m,r, phiVnm, phiWnm, phiXnm
// output: function values: fVnm,fWnm,fXnm,
//         function derivatives: fpVnm,fpWnm,fpXnm,
template <class Real, class Integer>
inline void fnExt(const Integer &n, const Integer &m, const Real &r, const Real &phiV, const Real &phiW,
                  const Real &phiX, Real &fVnm, Real &fWnm, Real &fXnm, Real &fpVnm, Real &fpWnm, Real &fpXnm) {
    fVnm = (Power(r, -2 - n) *
            ((2 * n * phiV) / (3 + 4 * n * (2 + n)) - ((1 + n) * phiW * (-1 + Power(r, 2))) / (1 + 2 * n))) /
           2.0;
    fWnm = ((1 + n) * phiW) / ((-1 + 4 * Power(n, 2)) * Power(r, n));
    fXnm = (phiX * Power(r, -1 - n)) / (1 + 2 * n);

    fpVnm = (Power(r, -3 - n) * ((-2 * n * (2 + n) * phiV) / (3 + 4 * n * (2 + n)) +
                                 ((1 + n) * phiW * (-2 + n * (-1 + Power(r, 2)))) / (1 + 2 * n))) /
            2.0;
    fpWnm = -((n * (1 + n) * phiW * Power(r, -1 - n)) / (-1 + 4 * Power(n, 2)));
    fpXnm = -(((1 + n) * phiX * Power(r, -2 - n)) / (1 + 2 * n));

    return;
}

template <class Real>
inline void fnInt(const Integer &n, const Integer &m, const Real &r, const Real &phiV, const Real &phiW,
                  const Real &phiX, Real &fVnm, Real &fWnm, Real &fXnm, Real &fpVnm, Real &fpWnm, Real &fpXnm) {

    fVnm = (n * phiV * Power(r, 1 + n)) / (3 + 4 * n * (2 + n));
    fWnm = -((1 + n) * Power(r, -1 + n) * (-2 * phiW + (-1 + 2 * n) * phiV * (-1 + Power(r, 2)))) /
           (2. * (-1 + 4 * Power(n, 2)));
    fXnm = (phiX * Power(r, n)) / (1 + 2 * n);

    fpVnm = (n * (1 + n) * phiV * Power(r, n)) / (3 + 4 * n * (2 + n));
    fpWnm = -((1 + n) * Power(r, -2 + n) *
              (-2 * (-1 + n) * phiW + (-1 + 2 * n) * phiV * (1 + Power(r, 2) + n * (-1 + Power(r, 2))))) /
            (2. * (-1 + 4 * Power(n, 2)));
    fpXnm = (n * phiX * Power(r, -1 + n)) / (1 + 2 * n);

    return;
}

#endif