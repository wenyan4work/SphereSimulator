/*
 *  Created on: Oct 6, 2016
 *      Author: wyan
 */

#ifndef INCLUDE_STKFMM_H
#define INCLUDE_STKFMM_H

#include <unordered_map>

// a wrapper for pvfmm
// choose kernel at compile time
#include <mpi.h>
#include <pvfmm.hpp>

class STKFMM {
  public:
    enum PAXIS {
        NONE,
        // TODO: add periodic BC
        PXYZ,
        PZ,
        PXY
    };

    enum class KERNEL : size_t {
        PVel = 1, // single layer kernel
        PVelGrad = 2,
        PVelLaplacian = 4,
        Traction = 8,
        LAPPGrad = 16, // laplace single layer
    };

    template <typename Enumeration>
    auto asInteger(Enumeration const value) -> typename std::underlying_type<Enumeration>::type {
        return static_cast<typename std::underlying_type<Enumeration>::type>(value);
    }

    STKFMM(int multOrder = 10, int maxPts = 1000, PAXIS pbc_ = PAXIS::NONE, unsigned int kernelComb_ = 1);

    ~STKFMM();

    void evaluateFMM(std::vector<double> &srcSLValue, std::vector<double> &srcDLValue, std::vector<double> &trgValue,
                     KERNEL kernelChoice);

    void clearFMM(KERNEL kernelChoice);

    void setPoints(const std::vector<double> &srcSLCoord_, const std::vector<double> &srcDLCoord_,
                   const std::vector<double> &trgCoord_);

    void setupTree(KERNEL kernel_);

    void setBox(double, double, double, double, double, double);

    void showActiveKernels();

    bool isKernelActive(KERNEL kernel_) { return asInteger(kernel_) & kernelComb; }

    void getKernelDimension(int &kdimSL_, int &kdimDL_, int &kdimTrg_, KERNEL kernel_) {
        kdimSL_ = poolFMM[kernel_]->kdimSL;
        kdimDL_ = poolFMM[kernel_]->kdimDL;
        kdimTrg_ = poolFMM[kernel_]->kdimTrg;
    }

  private:
    class FMMData {
      public:
        const KERNEL kernelChoice;
        const int multOrder;
        const size_t maxPts;

        int kdimSL;
        int kdimDL;
        int kdimTrg;

        // forbid default constructor
        FMMData() = delete;
        // copy constructor
        FMMData(const FMMData &) = delete;
        FMMData &operator=(const FMMData &) = delete;
        FMMData(FMMData &&) = delete;
        FMMData &operator=(FMMData &&) = delete;

        // constructor with choice of kernel
        FMMData(KERNEL, int, int);

        // destructor
        ~FMMData();

        // helper
        void setKernel(const pvfmm::Kernel<double> &kernelFunction);

        // computation routines
        void setupTree(const std::vector<double> &srcSLCoord, const std::vector<double> &srcDLCoord,
                       const std::vector<double> &trgCoord);
        void evaluate(std::vector<double> &srcSLValue, std::vector<double> &srcDLValue, std::vector<double> &trgValue);

        void deleteTree();
        void clear();

      private:
        pvfmm::PtFMM *matrixPtr;
        pvfmm::PtFMM_Tree *treePtr;
        pvfmm::PtFMM_Data *treeDataPtr;
    };

    const int multOrder;
    const int maxPts;
    PAXIS pbc;
    const unsigned int kernelComb;

    double xlow, xhigh; // box
    double ylow, yhigh;
    double zlow, zhigh;
    double scaleFactor;
    double xshift, yshift, zshift;

    std::vector<double> srcSLCoordInternal; // scaled coordinate Single Layer
    std::vector<double> srcDLCoordInternal; // scaled coordinate Double Layer
    std::vector<double> trgCoordInternal;
    std::vector<double> srcSLValueInternal; // scaled SL value
    std::vector<double> srcDLValueInternal; // scaled SL value
    std::vector<double> trgValueInternal;   // scaled trg value

    void setupCoord(const std::vector<double> &, std::vector<double> &); // setup the internal srcCoord and
                                                                         // trgCoord, with proper rotation and BC
    struct EnumClassHash {
        template <typename T>
        std::size_t operator()(T t) const {
            return static_cast<std::size_t>(t);
        }
    };

    std::unordered_map<KERNEL, STKFMM::FMMData *, EnumClassHash> poolFMM;
};

#endif /* INCLUDE_FMMWRAPPER_H_ */
