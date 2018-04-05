#ifndef GRADUNEAREVAL_HPP_
#define GRADUNEAREVAL_HPP_

#include "SLVel.hpp"

#include <complex>

template <class Integer>
inline bool validnm(const Integer &n, const Integer &m) {
    return (m <= n) && (m >= -n);
}

// input: n,m ; r,theta,phi; phiV, phiW, phiX
// output: 3x3 tensor gradVel, must be real
// TODO: check positive and negative m terms combined together in a single call ?
// g 0 1 2 are in r,theta,phi basis
template <class Real, class Integer>
void gradVelNearEval(const Integer &n, const Integer &m, const Real &r, const Real &theta, const Real &phi,
                     const Real &phiV, const Real &phiW, const Real &phiX, const bool &exterior, Real &g00, Real &g01,
                     Real &g02, Real &g10, Real &g11, Real &g12, Real &g20, Real &g21, Real &g22) {
    using Complex = std::complex<Real>;

    // TODO: replace the 1 with the actual value of Ynm
    const Complex Ynmtp = validnm(n, m) ? 1 : 0;
    const Complex Ynm1tp = validnm(n, m + 1) ? 1 : 0;
    const Complex Ynm2tp = validnm(n, m + 2) ? 1 : 0;

    Real fVnm = 0, fWnm = 0, fXnm = 0;
    Real fpVnm = 0, fpWnm = 0, fpXnm = 0;
    if (exterior) {
        fnExt(n, m, r, phiV, phiW, phiX, fVnm, fWnm, fXnm, fpVnm, fpWnm, fpXnm);
    } else {
        fnInt(n, m, r, phiV, phiW, phiX, fVnm, fWnm, fXnm, fpVnm, fpWnm, fpXnm);
    }

    Complex G00 = (fpVnm * (-1 - n) + fpWnm * n) * Ynmtp;
    Complex G01 =
        (-((Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * (fWnm - fWnm * n + fVnm * (2 + n)) * Ynm1tp) / Power(E, I * phi)) +
         fWnm * m * (-1 + n) * Ynmtp * Cot(theta) - fVnm * m * (2 + n) * Ynmtp * Cot(theta) +
         Complex(0, 1) * fXnm * m * Ynmtp * Csc(theta)) /
        r;
    Complex G02 = (-(fXnm * Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp) +
                   m * Ynmtp * (fWnm - fWnm * n + fVnm * (2 + n) - Complex(0, 1) * fXnm * Cos(theta)) * Csc(theta) *
                       (Complex(0, -1) * Cos(phi) + Sin(phi))) /
                  (Power(E, I * phi) * r);

    Complex G10 =
        (fpVnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp + fpWnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp +
         Power(E, I * phi) * m * Ynmtp * (Complex(0, -1) * fpXnm + (fpVnm + fpWnm) * Cos(theta)) * Csc(theta)) /
        Power(E, I * phi);
    Complex G11 =
        (Complex(0, -1) * Power(E, I * phi) * fXnm * m *
             (Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp + Power(E, I * phi) * (-1 + m) * Ynmtp * Cot(theta)) *
             Csc(theta) +
         fWnm * (Sqrt((m - n) * (1 + m - n) * (1 + m + n) * (2 + m + n)) * Ynm2tp +
                 Power(E, I * phi) * (1 + 2 * m) * Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp * Cot(theta) +
                 Power(E, 2 * I * phi) * Power(m, 2) * Ynmtp * Power(Cot(theta), 2) +
                 Power(E, 2 * I * phi) * Ynmtp * (n - m * Power(Csc(theta), 2))) +
         fVnm * (Sqrt((m - n) * (1 + m - n) * (1 + m + n) * (2 + m + n)) * Ynm2tp +
                 Power(E, I * phi) * (1 + 2 * m) * Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp * Cot(theta) +
                 Power(E, 2 * I * phi) * Power(m, 2) * Ynmtp * Power(Cot(theta), 2) -
                 Power(E, 2 * I * phi) * Ynmtp * (1 + n + m * Power(Csc(theta), 2)))) /
        (Power(E, 2 * I * phi) * r);
    Complex G12 =
        (-(fXnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Cot(theta)) +
         Complex(0, 1) * fVnm * m * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Csc(theta) +
         Complex(0, 1) * fWnm * m * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Csc(theta) +
         Power(E, I * phi) * m * Ynmtp *
             (fXnm + (-1 + m) * Csc(theta) * (Complex(0, 1) * (fVnm + fWnm) * Cot(theta) + fXnm * Csc(theta)))) /
        (Power(E, I * phi) * r);

    Complex G20 = (fpXnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp) / Power(E, I * phi) +
                  m * Ynmtp * (Complex(0, 1) * (fpVnm + fpWnm) + fpXnm * Cos(theta)) * Csc(theta);
    Complex G21 =
        (Complex(0, 1) * Power(E, I * phi) * fVnm * m *
             (Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp + Power(E, I * phi) * (-1 + m) * Ynmtp * Cot(theta)) *
             Csc(theta) +
         Complex(0, 1) * Power(E, I * phi) * fWnm * m *
             (Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp + Power(E, I * phi) * (-1 + m) * Ynmtp * Cot(theta)) *
             Csc(theta) +
         (fXnm * (2 * Sqrt((m - n) * (1 + m - n) * (1 + m + n) * (2 + m + n)) * Ynm2tp +
                  2 * Power(E, I * phi) * (1 + 2 * m) * Sqrt(-(m * (1 + m)) + n + Power(n, 2)) * Ynm1tp * Cot(theta) +
                  Power(E, 2 * I * phi) * m * Ynmtp * (-2 + m + m * Cos(2 * theta)) * Power(Csc(theta), 2))) /
             2.) /
        (Power(E, 2 * I * phi) * r);
    Complex G22 = (fVnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Cot(theta) +
                   fWnm * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Cot(theta) +
                   Complex(0, 1) * fXnm * m * Sqrt((-m + n) * (1 + m + n)) * Ynm1tp * Csc(theta) -
                   Power(E, I * phi) * Ynmtp *
                       (Complex(0, -1) * fXnm * (-1 + m) * m * Cot(theta) * Csc(theta) +
                        fVnm * (1 + m + n + (-1 + m) * m * Power(Csc(theta), 2)) -
                        fWnm * (n + m * (Power(Cot(theta), 2) - m * Power(Csc(theta), 2))))) /
                  (Power(E, I * phi) * r);

    // TODO: find the relation between complex G and real g

    return;
}

#endif