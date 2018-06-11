/*
 * GeoUtil.h
 *
 *  Created on: Sep 19, 2016
 *      Author: wyan
 */

#ifndef TUBULESIMULATOR_TUBULESIMULATOR_GEOUTIL_H_
#define TUBULESIMULATOR_TUBULESIMULATOR_GEOUTIL_H_

//
// Created by Wen Yan on 9/19/16.
//

#include <cmath>
#include <limits>

template <typename Vector>
struct ResultP3S3 {
    double distance, sqrDistance;
    double segmentParameter; // t in [0,1]
    Vector segmentClosest;   // (1-t)*p[0] + t*p[1]
};

template <typename Vector>
struct ResultS3S3 {
    double distance, sqrDistance;
    double segmentParameter[2]; // t in [0,1]
    Vector segmentClosest[2];   // (1-t)*p[0] + t*p[1]
};

template <typename Vec3>
inline double norm(Vec3 &a) {
    return a.norm();
}

template <typename Vec3>
inline void normalize(Vec3 &a) {
    // double r2 = a * a;
    // a = a * (1 / sqrt(r2));
    a.normalize();
}

template <typename Vector>
inline double dot(Vector &a, Vector &b) {
    return a.dot(b);
}

template <typename Vector>
double DistPointSeg(const Vector &point, const Vector &minus, const Vector &plus, Vector &pointPerp) {
    ResultP3S3<Vector> result;

    // The direction vector is not unit length.  The normalization is deferred
    // until it is needed.
    Vector direction = plus - minus;
    Vector diff = point - plus;
    double t = dot(direction, diff);
    if (t >= 0) {
        result.segmentParameter = 1;
        result.segmentClosest = plus;
    } else {
        diff = point - minus;
        t = dot(direction, diff);
        if (t <= 0) {
            result.segmentParameter = 0;
            result.segmentClosest = minus;
        } else {
            double sqrLength = dot(direction, direction);
            if (sqrLength > 0) {
                t /= sqrLength;
                result.segmentParameter = t;
                result.segmentClosest = minus + t * direction;
            } else {
                result.segmentParameter = 0;
                result.segmentClosest = minus;
            }
        }
    }

    diff = point - result.segmentClosest;
    result.sqrDistance = dot(diff, diff);
    result.distance = sqrt(result.sqrDistance);

    pointPerp = result.segmentClosest;

    return result.distance;
}

// this object must be thread private in multi-threading environment
template <int N, typename Real, typename Vector>
class DCPQuery {
  public:
    struct Result {
        Real distance, sqrDistance;
        Real parameter[2];
        Vector closest[2];
    };

    Real operator()(Vector const &P0, Vector const &P1, Vector const &Q0, Vector const &Q1, Vector &Ploc, Vector &Qloc,
                    Real &s, Real &t);

    Real operator()(Vector const &P0, Vector const &P1, Vector const &Q0, Vector const &Q1, Vector &Ploc, Vector &Qloc);

  private:
    // Compute the root of h(z) = h0 + slope*z and clamp it to the interval
    // [0,1].  It is required that for h1 = h(1), either (h0 < 0 and h1 > 0)
    // or (h0 > 0 and h1 < 0).
    Real GetClampedRoot(Real slope, Real h0, Real h1);

    // Compute the intersection of the line dR/ds = 0 with the domain [0,1]^2.
    // The direction of the line dR/ds is conjugate to (1,0), so the algorithm
    // for minimization is effectively the conjugate gradient algorithm for a
    // quadratic function.
    void ComputeIntersection(Real const sValue[2], int const classify[2], int edge[2], Real end[2][2]);

    // Compute the location of the minimum of R on the segment of intersection
    // for the line dR/ds = 0 and the domain [0,1]^2.
    void ComputeMinimumParameters(int const edge[2], Real const end[2][2], Real parameter[2]);

    // The coefficients of R(s,t), not including the constant term.
    Real mA, mB, mC, mD, mE;

    // dR/ds(i,j) at the four corners of the domain
    Real mF00, mF10, mF01, mF11;

    // dR/dt(i,j) at the four corners of the domain
    Real mG00, mG10, mG01, mG11;
};

template <int N, typename Real, typename Vector>
Real DCPQuery<N, Real, Vector>::operator()(Vector const &P0, Vector const &P1, Vector const &Q0, Vector const &Q1,
                                           Vector &Ploc, Vector &Qloc) {
    Real s, t = 0;
    return (*this)(P0, P1, Q0, Q1, Ploc, Qloc, s, t);
}

template <int N, typename Real, typename Vector>
Real DCPQuery<N, Real, Vector>::operator()(Vector const &P0, Vector const &P1, Vector const &Q0, Vector const &Q1,
                                           Vector &Ploc, Vector &Qloc, Real &s, Real &t) {
    DCPQuery<N, Real, Vector>::Result result;

    // The code allows degenerate line segments; that is, P0 and P1 can be
    // the same point or Q0 and Q1 can be the same point.  The quadratic
    // function for squared distance between the segment is
    //   R(s,t) = a*s^2 - 2*b*s*t + c*t^2 + 2*d*s - 2*e*t + f
    // for (s,t) in [0,1]^2 where
    //   a = Dot(P1-P0,P1-P0), b = Dot(P1-P0,Q1-Q0), c = Dot(Q1-Q0,Q1-Q0),
    //   d = Dot(P1-P0,P0-Q0), e = Dot(Q1-Q0,P0-Q0), f = Dot(P0-Q0,P0-Q0)
    Vector P1mP0 = P1 - P0;
    Vector Q1mQ0 = Q1 - Q0;
    Vector P0mQ0 = P0 - Q0;
    mA = dot(P1mP0, P1mP0);
    mB = dot(P1mP0, Q1mQ0);
    mC = dot(Q1mQ0, Q1mQ0);
    mD = dot(P1mP0, P0mQ0);
    mE = dot(Q1mQ0, P0mQ0);

    mF00 = mD;
    mF10 = mF00 + mA;
    mF01 = mF00 - mB;
    mF11 = mF10 - mB;

    mG00 = -mE;
    mG10 = mG00 - mB;
    mG01 = mG00 + mC;
    mG11 = mG10 + mC;

    if (mA > (Real)0 && mC > (Real)0) {
        // Compute the solutions to dR/ds(s0,0) = 0 and dR/ds(s1,1) = 0.  The
        // location of sI on the s-axis is stored in classifyI (I = 0 or 1).  If
        // sI <= 0, classifyI is -1.  If sI >= 1, classifyI is 1.  If 0 < sI < 1,
        // classifyI is 0.  This information helps determine where to search for
        // the minimum point (s,t).  The fij values are dR/ds(i,j) for i and j in
        // {0,1}.

        Real sValue[2];
        sValue[0] = GetClampedRoot(mA, mF00, mF10);
        sValue[1] = GetClampedRoot(mA, mF01, mF11);

        int classify[2];
        for (int i = 0; i < 2; ++i) {
            if (sValue[i] <= (Real)0) {
                classify[i] = -1;
            } else if (sValue[i] >= (Real)1) {
                classify[i] = +1;
            } else {
                classify[i] = 0;
            }
        }

        if (classify[0] == -1 && classify[1] == -1) {
            // The minimum must occur on s = 0 for 0 <= t <= 1.
            result.parameter[0] = (Real)0;
            result.parameter[1] = GetClampedRoot(mC, mG00, mG01);
        } else if (classify[0] == +1 && classify[1] == +1) {
            // The minimum must occur on s = 1 for 0 <= t <= 1.
            result.parameter[0] = (Real)1;
            result.parameter[1] = GetClampedRoot(mC, mG10, mG11);
        } else {
            // The line dR/ds = 0 intersects the domain [0,1]^2 in a
            // nondegenerate segment.  Compute the endpoints of that segment,
            // end[0] and end[1].  The edge[i] flag tells you on which domain
            // edge end[i] lives: 0 (s=0), 1 (s=1), 2 (t=0), 3 (t=1).
            int edge[2];
            Real end[2][2];
            ComputeIntersection(sValue, classify, edge, end);

            // The directional derivative of R along the segment of
            // intersection is
            //   H(z) = (end[1][1]-end[1][0])*dR/dt((1-z)*end[0] + z*end[1])
            // for z in [0,1].  The formula uses the fact that dR/ds = 0 on
            // the segment.  Compute the minimum of H on [0,1].
            ComputeMinimumParameters(edge, end, result.parameter);
        }
    } else {
        if (mA > (Real)0) {
            // The Q-segment is degenerate (Q0 and Q1 are the same point) and
            // the quadratic is R(s,0) = a*s^2 + 2*d*s + f and has (half)
            // first derivative F(t) = a*s + d.  The closest P-point is
            // interior to the P-segment when F(0) < 0 and F(1) > 0.
            result.parameter[0] = GetClampedRoot(mA, mF00, mF10);
            result.parameter[1] = (Real)0;
        } else if (mC > (Real)0) {
            // The P-segment is degenerate (P0 and P1 are the same point) and
            // the quadratic is R(0,t) = c*t^2 - 2*e*t + f and has (half)
            // first derivative G(t) = c*t - e.  The closest Q-point is
            // interior to the Q-segment when G(0) < 0 and G(1) > 0.
            result.parameter[0] = (Real)0;
            result.parameter[1] = GetClampedRoot(mC, mG00, mG01);
        } else {
            // P-segment and Q-segment are degenerate.
            result.parameter[0] = (Real)0;
            result.parameter[1] = (Real)0;
        }
    }

    result.closest[0] = ((Real)1 - result.parameter[0]) * P0 + result.parameter[0] * P1;
    result.closest[1] = ((Real)1 - result.parameter[1]) * Q0 + result.parameter[1] * Q1;
    Vector diff = result.closest[0] - result.closest[1];
    result.sqrDistance = dot(diff, diff);
    result.distance = sqrt(result.sqrDistance);
    Ploc = result.closest[0];
    Qloc = result.closest[1];
    s = result.parameter[0];
    t = result.parameter[1];
    return result.distance;
}

template <int N, typename Real, typename Vector>
inline Real DCPQuery<N, Real, Vector>::GetClampedRoot(Real slope, Real h0, Real h1) {
    // slope = h1-h0
    // h0 and h1 should have different sign, return the zero point
    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    // if (slope - (h1 - h0) > eps) {
    //     std::cout << "slope" << slope << std::endl;
    //     std::cout << "h1" << h1 << std::endl;
    //     std::cout << "h0" << h0 << std::endl;
    //     std::cout << "diff" << slope - (h1 - h0) << std::endl;
    // }
    assert(fabs(slope - (h1 - h0)) < 10 * eps * std::max(abs(h1), abs(h0)));

    Real r;

    if (std::abs(h0) < eps && std::abs(h1) < eps) {
        // tiny slope, h0 \approx h1, distance almost a constant, choose mid point
        r = 0.5;
    } else if (h0 < 0) {
        if (h1 > 0) {
            // r = -h0 / slope; // need better accuracy
            // clamp r between [0,1]
            r = std::min(std::max(-h0 / slope, (Real)0), (Real)1);
        } else {
            r = 1;
        }
    } else {
        r = 0;
    }
    return r;
}

template <int N, typename Real, typename Vector>
inline void DCPQuery<N, Real, Vector>::ComputeIntersection(Real const sValue[2], int const classify[2], int edge[2],
                                                           Real end[2][2]) {
    // The divisions are theoretically numbers in [0,1].  Numerical rounding
    // errors might cause the result to be outside the interval.  When this
    // happens, it must be that both numerator and denominator are nearly
    // zero.  The denominator is nearly zero when the segments are nearly
    // perpendicular.  The numerator is nearly zero when the P-segment is
    // nearly degenerate (mF00 = a is small).  The choice of 0.5 should not
    // cause significant accuracy problems.
    //
    // NOTE:  You can use bisection to recompute the root or even use
    // bisection to compute the root and skip the division.  This is generally
    // slower, which might be a problem for high-performance applications.

    if (classify[0] < 0) {
        edge[0] = 0;
        end[0][0] = (Real)0;
        end[0][1] = mF00 / mB;
        if (end[0][1] < (Real)0 || end[0][1] > (Real)1) {
            end[0][1] = (Real)0.5;
        }

        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = (Real)1;
        } else // classify[1] > 0
        {
            edge[1] = 1;
            end[1][0] = (Real)1;
            end[1][1] = mF10 / mB;
            if (end[1][1] < (Real)0 || end[1][1] > (Real)1) {
                end[1][1] = (Real)0.5;
            }
        }
    } else if (classify[0] == 0) {
        edge[0] = 2;
        end[0][0] = sValue[0];
        end[0][1] = (Real)0;

        if (classify[1] < 0) {
            edge[1] = 0;
            end[1][0] = (Real)0;
            end[1][1] = mF00 / mB;
            if (end[1][1] < (Real)0 || end[1][1] > (Real)1) {
                end[1][1] = (Real)0.5;
            }
        } else if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = (Real)1;
        } else {
            edge[1] = 1;
            end[1][0] = (Real)1;
            end[1][1] = mF10 / mB;
            if (end[1][1] < (Real)0 || end[1][1] > (Real)1) {
                end[1][1] = (Real)0.5;
            }
        }
    } else // classify[0] > 0
    {
        edge[0] = 1;
        end[0][0] = (Real)1;
        end[0][1] = mF10 / mB;
        if (end[0][1] < (Real)0 || end[0][1] > (Real)1) {
            end[0][1] = (Real)0.5;
        }

        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = (Real)1;
        } else {
            edge[1] = 0;
            end[1][0] = (Real)0;
            end[1][1] = mF00 / mB;
            if (end[1][1] < (Real)0 || end[1][1] > (Real)1) {
                end[1][1] = (Real)0.5;
            }
        }
    }
}

template <int N, typename Real, typename Vector>
inline void DCPQuery<N, Real, Vector>::ComputeMinimumParameters(int const edge[2], Real const end[2][2],
                                                                Real parameter[2]) {
    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    Real delta = end[1][1] - end[0][1];
    Real h0 = delta * ((-mB * end[0][0] - mE) + mC * end[0][1]); // source of rounding error
    Real h1 = delta * ((-mB * end[1][0] - mE) + mC * end[1][1]);

    if (std::abs(h0) < std::abs(mC) * eps && std::abs(h1) < std::abs(mC) * eps) {
        Real z = 0.5;
        Real omz = (Real)1 - z;
        parameter[0] = omz * end[0][0] + z * end[1][0];
        parameter[1] = omz * end[0][1] + z * end[1][1];
    } else if (h0 >= (Real)0) {
        if (edge[0] == 0) {
            parameter[0] = (Real)0;
            parameter[1] = GetClampedRoot(mC, mG00, mG01);
        } else if (edge[0] == 1) {
            parameter[0] = (Real)1;
            parameter[1] = GetClampedRoot(mC, mG10, mG11);
        } else {
            parameter[0] = end[0][0];
            parameter[1] = end[0][1];
        }
    } else {
        if (h1 <= (Real)0) {
            if (edge[1] == 0) {
                parameter[0] = (Real)0;
                parameter[1] = GetClampedRoot(mC, mG00, mG01);
            } else if (edge[1] == 1) {
                parameter[0] = (Real)1;
                parameter[1] = GetClampedRoot(mC, mG10, mG11);
            } else {
                parameter[0] = end[1][0];
                parameter[1] = end[1][1];
            }
        } else // h0 < 0 and h1 > 0
        {
            Real z = GetClampedRoot(h1 - h0, h0, h1);
            Real omz = (Real)1 - z;
            parameter[0] = omz * end[0][0] + z * end[1][0];
            parameter[1] = omz * end[0][1] + z * end[1][1];
        }
    }

    return;
}

// Template aliases for convenience.
template <int N, typename Real, typename Vector>
using DCPSegmentSegment = DCPQuery<N, Real, Vector>;

template <typename Real, typename Vector>
using DCPSegment2Segment2 = DCPSegmentSegment<2, Real, Vector>;

template <typename Real, typename Vector>
using DCPSegment3Segment3 = DCPSegmentSegment<3, Real, Vector>;

//
///*Geometric Utility functions for use*/
// template<typename T >
// double DistPointSeg(const T point, const T segMinus, const T segPlus,
//		PS::F64vec3 & pointPerp) {
//	const auto vec1 = (point - segMinus) ^ (point - segPlus);
//	const auto vec2 = (segPlus - segMinus);
//	// use orthogonal projection to get the pointPerp
//	const auto vec3 = (point - segMinus);
//	pointPerp = segMinus + (vec3 * vec2) * (vec2 / (vec2 * vec2));
//
//	return sqrt(vec1 * vec1 / (vec2 * vec2));
//}

// the original function
// template <int N, typename Real, typename Vector>
// inline Real DCPQuery<N, Real, Vector>::GetClampedRoot(Real slope, Real h0, Real h1) {
//     // Theoretically, r is in (0,1).  However, when the slope is nearly zero,
//     // then so are h0 and h1.  Significant numerical rounding problems can
//     // occur when using floating-point arithmetic.  If the rounding causes r
//     // to be outside the interval, clamp it.  It is possible that r is in
//     // (0,1) and has rounding errors, but because h0 and h1 are both nearly
//     // zero, the quadratic is nearly constant on (0,1).  Any choice of p
//     // should not cause undesirable accuracy problems for the final distance
//     // computation.
//     //
//     // NOTE:  You can use bisection to recompute the root or even use
//     // bisection to compute the root and skip the division.  This is generally
//     // slower, which might be a problem for high-performance applications.

//     Real r;

//  if (h0 < (Real)0) {
//         if (h1 > (Real)0) {
//             r = -h0 / slope;
//             if (r > (Real)1) {
//                 r = (Real)0.5;
//             }
//             // The slope is positive and -h0 is positive, so there is no
//             // need to test for a negative value and clamp it.
//         } else {
//             r = (Real)1;
//         }
//     } else {
//         r = (Real)0;
//     }
//     return r;
// }

#endif /* TUBULESIMULATOR_TUBULESIMULATOR_GEOUTIL_H_ */
