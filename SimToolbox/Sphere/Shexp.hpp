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

// naming convention:
// Spectral Coeff
// Grid Value

// data arrangement:
// Spectral Coeff SHCArrange::Row_Major, coeff for X component, coeff for Y component, coeff for Z component
// Grid Value (x0,y0,z0,x1,y1,z1,x2,y2,z2,....)

// scaling of radius:
// No scaling of radius in calcGridValue / calcPoleValue / calcSpectralCoeff functions
// Correct scaling in near eval routines

class Shexp {
  public:
    enum class KIND {
        LAP, // Laplace single layer, etc
        STK  // Stokes
    };

    KIND kind;
    int order;          // the order of expansion p
    std::string name;   // the name of this quantity
    Equatn orientation; // the orientation represented by an Eigen::quaternion. Z axis = orientation * Evec3(0,0,1)
    double radius;      // radius of this Sph

    std::vector<double> gridValue;
    // dimension real numbers per point, (p+1) * (2p+2) points, excluding the north and south poles

    // constructor
    Shexp(const KIND kind_, const std::string &name_, const int order_, const double &radius = 1,
          const Equatn orientation_ = Equatn::Identity());

    // copy constructor
    Shexp(const Shexp &) = default;
    Shexp &operator=(const Shexp &) = default;

    // move constructor
    Shexp(Shexp &&);
    Shexp &operator=(Shexp &&);

    // destructor
    ~Shexp() = default;

    // utility routines
    inline int getDimension() const { return kind == KIND::LAP ? 1 : 3; }
    inline int getGridDOF() const { return getDimension() * getGridNumber(); }
    inline int getGridNumber() const { return (order + 1) * (2 * order + 2); }
    inline int getSpectralDOF() const { return getDimension() * (order + 1) * (order + 2); }

    // grid representation
    void getGridWithPole(std::vector<double> &gridPoints, std::vector<double> &gridWeights,
                         const Evec3 &coordBase = Evec3::Zero(), std::vector<double> *gridNormPtr = nullptr) const;

    // output routine
    // each sph dump to a 'piece'.
    // points in different pieces are completely independent
    // variable names in different pieces must be the same.
    int writeVTU(std::ofstream &file, const Evec3 &coordBase = Evec3::Zero()) const;

    // debug routines
    void dumpSpectralCoeff(const double *const spectralCoeff, const std::string &filename = std::string("")) const;
    void dumpGridValue(const std::string &filename = std::string("")) const;
    void randomFill(const int seed = 0);

    // convert G2S and S2G
    // User should allocate enough space. No bound check here
    void calcGridValue(double *coeffPtr, double *valPtr = nullptr);
    void calcPoleValue(double *coeffPtr, double *valPtr) const;

    void calcSpectralCoeff(double *coeffPtr, double *valPtr = nullptr) const;

    // NF routines
    // User should allocate enough space. No bound check here
    void calcSDLNF(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgValuePtr,
                   const bool &interior = false, const bool &SL = true) const; // both STK and LAP

    void calcKNF(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgNormPtr, double *trgValuePtr,
                 const bool &interior = false) const;

    void calcKSelf(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgValuePtr,
                   const bool &interior = false) const;
    // points are all on a sphere centered at the sph grid center, norm parallel with XYZ, always outward

  private:
    // called from dumpVTK()
    void getGridWithPoleCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
                                    std::vector<uint8_t> &type) const;

    void rotGridValue(double *valPtr, const int npts) const;    // from default Z axis to oriented Z axis
    void invrotGridValue(double *valPtr, const int npts) const; // from oriented Z axis to default Z axis
};

#endif // SPHEXP_HPP
