/*
 * FMMWrapper.h
 *
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
        SLPVel = 1, // single layer kernel
        SLPVelGrad = 2,
        SLTraction = 4,
        SLPVelLaplacian = 8,
        DLPVel = 16, // double layer kernel
        DLPVelGrad = 32,
        LAPSLPGrad = 64, // laplace single layer
        LAPDLPGrad = 128, // laplace double layer
    };

    template <typename Enumeration>
    auto asInteger(Enumeration const value) -> typename std::underlying_type<Enumeration>::type {
        return static_cast<typename std::underlying_type<Enumeration>::type>(value);
    }

    STKFMM(int multOrder = 10, int maxPts = 1000, PAXIS pbc_ = PAXIS::NONE, unsigned int kernelComb_ = 1);

    ~STKFMM();

    void evaluateFMM(std::vector<double> &srcValue, std::vector<double> &trgValue, KERNEL kernelChoice);

    void clearFMM(KERNEL kernelChoice);

    void setPoints(const std::vector<double> &srcCoord_, const std::vector<double> &trgCoord_);

    void setupTree(KERNEL);

    void setBox(double, double, double, double, double, double);

    void showActiveKernels();

    bool isKernelActive(KERNEL kernel_) { return asInteger(kernel_) & kernelComb; }

    void getKernelDimension(int &kdimSrc_, int &kdimTrg_, KERNEL kernel_) {
        kdimSrc_ = poolFMM[kernel_]->kdimSrc;
        kdimTrg_ = poolFMM[kernel_]->kdimTrg;
    }

  private:
    class FMMData {
      public:
        const KERNEL kernelChoice;
        const int multOrder;
        const size_t maxPts;

        int kdimSrc;
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
        void setupTree(const std::vector<double> &, const std::vector<double> &);
        void deleteTree();
        void clear();
        void evaluate(std::vector<double> &, std::vector<double> &);

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

    std::vector<double> srcCoordInternal; // scaled coordinate
    std::vector<double> trgCoordInternal;
    // std::vector<double> srcValueInternal;
    std::vector<double> trgValueInternal; // scaled trg value

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
