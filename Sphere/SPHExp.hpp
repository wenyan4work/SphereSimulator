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
    int order;          // the order of expansion p
    std::string name;   // the name of this quantity
    Equatn orientation; // the orientation represented by an Eigen::quaternion

    std::vector<double>
        gridValues; // dimension real numbers per point, (p+1) * (2p+2)+2 points, including the north and south pole

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
    inline int getGridDOF() const { return getDimension() * (order + 1) * (2 * order + 2); }
    inline int getSpectralDOF() const { return getDimension() * (order + 1) * (order + 2); }

    // grid representation
    void getGrid(std::vector<double> &gridPoints, std::vector<double> &gridWeights, const double &radius = 1,
                 const Evec3 &coordBase = Evec3::Zero()) const;

    // output routine
    // each sph dump to a 'piece'.
    // points in different pieces are completely independent
    // variable names in different pieces must be the same.
    int writeVTU(std::ofstream &file, const double &radius = 1, const Evec3 &coordBase = Evec3::Zero()) const;

    // debug routines
    void dumpSpectralValues(const double *constspectralCoeff,
                            const std::string &filename = std::string("")) const; // default to empty string

    // convert G2S and S2G
    // User should allocate enough space. No bound check here
    void calcGridValue(double *valPtr, double *coeffPtr) const;
    void calcGridValuePole(double *valPtr, double *coeffPtr) const;

    void calcSpectralValue(double *coeffPtr, double *valPtr) const;

    // NF routines
    // User should allocate enough space. No bound check here
    void calcSLNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const;   // both STK and LAP
    void calcDLNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const;   // both STK and LAP
    void calcTracNF(double *trgValue, const double *const trgGrid, const TRGPOS &trgPos) const; // STK Traction

  private:
    // called from dumpVTK()
    void getGridCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
                            std::vector<uint8_t> &type) const;
};

#endif // SPHEXP_HPP
