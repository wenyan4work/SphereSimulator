#ifndef SPHEXP_HPP
#define SPHEXP_HPP

// the class for spherical harmonics expansions
// should be guaranteed to be thread-safe
#include <array>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "Util/EigenDef.hpp"

class SPHExp {
  public:
    enum class KIND {
        LAP, // Laplace single layer, etc
        STK  // Stokes
    };

    enum class TRGPOS {
        IN, // target inside, outside, on the sphere
        OUT
    };

    KIND kind;
    int order; // the order of expansion p
    // int dimension;      // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL
    std::string name;   // the name of this quantity
    Equatn orientation; // the orientation represented by an Eigen::quaternion

    std::vector<double> spectralCoeff; // coefficients. d* p*(p+1) elements

    // index
    inline int COEFFINDEX(int i, int n, int m) const {
        // i=0 for LAPSL and LAPDL
        // i=0,1,2 for STKSL and STKDL
        return i * (n + 1) * (n + 1) + n * n + m + n;
    }

    // constructor
    SPHExp(const KIND kind_, const std::string &name_, const int order_,
           const Equatn orientation_ = Equatn::Identity());

    // copy constructor
    SPHExp(const SPHExp &) = default;
    SPHExp &operator=(const SPHExp &) = default;

    // move constructor
    SPHExp(SPHExp &&);
    SPHExp &operator=(SPHExp &&);

    // destructor
    ~SPHExp() = default;

    // utility routines
    inline int getDimension() const { return kind == KIND::LAP ? 1 : 3; }

    // grid representation
    void getGrid(std::vector<double> &gridPoints, std::vector<double> &gridWeights, std::vector<double> &gridValues,
                 const double &radius = 1, const Evec3 &coordBase = Evec3::Zero()) const;

    int getGridDOF() const; // excluding the north and south pole
    int getSpectralDOF() const;

    // output routine
    // each sph dump to a 'piece'.
    // points in different pieces are completely independent
    // variable names in different pieces must be the same.
    int writeVTU(std::ofstream &file, const double &radius = 1, const Evec3 &coordBase = Evec3::Zero()) const;

    // debug routines
    void dumpSpectralValues(const std::string &filename = std::string("")) const; // default to empty string

    void getGridCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
                            std::vector<uint8_t> &type) const;

    // convert G2S and S2G
    void calcGridValues(double *val, const double *const coeff) const;
    void calcSpectralValues(double *coeff, const double *const val) const;

    // NF routines
    void calcSLNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const;   // both STK and LAP
    void calcDLNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const;   // both STK and LAP
    void calcTracNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const; // STK Traction
};

#endif // SPHEXP_HPP
