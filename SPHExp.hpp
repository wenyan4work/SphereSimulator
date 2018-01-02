#ifndef SPHEXP_HPP
#define SPHEXP_HPP

// the class for spherical harmonics expansions
#include <string>
#include <vector>
#include "EigenDef.hpp"

class SPHExp {
  public:
    enum KIND {
        LAPSL, // laplace single layer, etc
        LAPDL,
        STKSL,
        STKDL
    };

    enum TRGPOS {
        IN, // target inside, outside, on the sphere
        OUT,
        ON,
        CHECK // check with given r
    };

    const KIND kind;
    const int order;        // the order of expansion p
    const int dimension;    // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL
    const std::string name; // the name of this quantity

    std::vector<double> spectralCoeff; // coefficients. d* p*(p+1) elements
    Equatn orientation;                // the orientation represented by an Eigen::quaternion

    // index
    inline int COEFFINDEX(int i, int n, int m) {
        // i=0 for LAPSL and LAPDL
        // i=0,1,2 for STKSL and STKDL
        return i * (n + 1) * (n + 1) + n * n + m + n;
    }

    // constructor
    SPHExp(const KIND kind_, const std::string &name_, const int order_,
           const Equatn orientation_ = Equatn::Identity());

    // copy constructor
    SPHExp(const SPHExp &);
    SPHExp(SPHExp &&);

    // copy assignment
    SPHExp &operator=(SPHExp &) = delete;
    SPHExp &operator=(SPHExp &&) = delete;

    // destructor
    ~SPHExp() = default;

    // utility routines
    //    inline double &operator[](int k) { return spectralCoeff[k]; } // return the located at the given index
    //    inline double &operator()(int i, int n, int m) { return spectralCoeff[COEFFINDEX(i, n, m)]; }

    // orientation
    void rotateOrientation(const Evec3 &);
    void setOrientation(const Equatn &);

    // grid representation
    void getGridPoints(std::vector<double> &) const;
    void getGridValues(std::vector<double> &) const;
    void calcGridValues(std::vector<double> &, const std::vector<double> &) const;

    // Spectral representation
    void calcSpectralValues(std::vector<double> &, const std::vector<double> &) const;
    void setSpectralValues(const std::vector<double> &);
    void calcAndSetSpectralValues(const std::vector<double> &);

    // FF routines
    void evaluateSphCoord(std::vector<double> &, const std::vector<double> &, TRGPOS trgPos);
    void evaluateCartCoord(std::vector<double> &, const std::vector<double> &, TRGPOS trgPos);
    void integralSphCoord(std::vector<double> &, const std::vector<double> &, TRGPOS trgPos);
    void integralCartCoord(std::vector<double> &, const std::vector<double> &, TRGPOS trgPos);

    // NF  routines

    // debug routines
    void dumpVTK(const std::string &);
    void dumpSpectralValues(const std::string &filename = std::string("")); // default to empty string
};

#endif // SPHEXP_HPP
