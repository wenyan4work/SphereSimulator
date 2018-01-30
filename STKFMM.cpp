/*
 * FMMWrapper.cpp
 *
 *  Created on: Oct 20, 2016
 *      Author: wyan
 */

#include <bitset>
#include <cassert>

#include <mpi.h>
#include <omp.h>

#include "../common/Timer.hpp"
#include "LaplaceLayerKernel.hpp"
#include "StokesLayerKernel.hpp"

#include "STKFMM.h"

PeriodicType periodicType;

// TreeClear(), optional, should be called BEFORE TreeUpdate.

// return fraction part between [0,1)
/*
 * This function is only applied in the PERIODIC DIRECTION
 * The user of the library must ensure that all points are located within [0,1)
 * */
inline double fracwrap(double x) { return x - floor(x); }

template <class T>
void safeDeletePtr(T *ptr) {
    if (ptr != nullptr) {
        delete ptr;
        ptr = nullptr;
    }
}

void STKFMM::FMMData::setKernel(const pvfmm::Kernel<double> &kernelFunction) {
    matrixPtr->Initialize(multOrder, MPI_COMM_WORLD, &kernelFunction);
    kdimSL = kernelFunction.k_s2t->ker_dim[0];
    kdimTrg = kernelFunction.k_s2t->ker_dim[1];
    kdimDL = kernelFunction.surf_dim;
}

// constructor
STKFMM::FMMData::FMMData(KERNEL kernelChoice_, int multOrder_, int maxPts_)
    : kernelChoice(kernelChoice_), multOrder(multOrder_), maxPts(maxPts_), treePtr(nullptr), matrixPtr(nullptr),
      treeDataPtr(nullptr) {
    matrixPtr = new pvfmm::PtFMM();
    // choose a kernel
    switch (kernelChoice) {
    case KERNEL::PVel:
        setKernel(pvfmm::StokesLayerKernel<double>::PVel());
        break;
    case KERNEL::PVelGrad:
        setKernel(pvfmm::StokesLayerKernel<double>::PVelGrad());
        break;
    case KERNEL::PVelLaplacian:
        setKernel(pvfmm::StokesLayerKernel<double>::PVelLaplacian());
        break;
    case KERNEL::Traction:
        setKernel(pvfmm::StokesLayerKernel<double>::Traction());
        break;
    case KERNEL::LAPPGrad:
        setKernel(pvfmm::LaplaceLayerKernel<double>::PGrad());
        break;
    }
    treeDataPtr = new pvfmm::PtFMM_Data;
    // treeDataPtr remain nullptr after constructor
}

STKFMM::FMMData::~FMMData() {
    clear();
    safeDeletePtr(treePtr);
    safeDeletePtr(treeDataPtr);
    safeDeletePtr(matrixPtr);
}

void STKFMM::FMMData::clear() {
    //    treeDataPtr->Clear();
    if (treePtr != nullptr)
        treePtr->ClearFMMData();
    return;
}

void STKFMM::FMMData::setupTree(const std::vector<double> &srcSLCoord, const std::vector<double> &srcDLCoord,
                                const std::vector<double> &trgCoord) {
    // trgCoord and srcCoord have been scaled to [0,1)^3

    // setup treeData
    treeDataPtr->dim = 3;
    treeDataPtr->max_depth = 15; // must < MAX_DEPTH in pvfmm_common.hpp
    treeDataPtr->max_pts = maxPts;

    treeDataPtr->src_coord = srcSLCoord;
    treeDataPtr->surf_coord = srcDLCoord;
    treeDataPtr->trg_coord = trgCoord;

    // this is used to setup FMM octree
    treeDataPtr->pt_coord = srcSLCoord.size() > trgCoord.size() ? srcSLCoord : trgCoord;
    const size_t nSL = srcSLCoord.size() / 3;
    const size_t nDL = srcDLCoord.size() / 3;
    const size_t nTrg = trgCoord.size() / 3;

    // space allocate
    treeDataPtr->src_value.Resize(nSL * kdimSL);
    treeDataPtr->surf_value.Resize(nDL * kdimDL);
    treeDataPtr->trg_value.Resize(nTrg * kdimTrg);

    // construct tree
    treePtr = new pvfmm::PtFMM_Tree(MPI_COMM_WORLD);
    treePtr->Initialize(treeDataPtr);
    treePtr->InitFMM_Tree(true, periodicType == PeriodicType::NONE ? pvfmm::FreeSpace : pvfmm::Periodic);
    treePtr->SetupFMM(matrixPtr);
    return;
}

void STKFMM::FMMData::deleteTree() {
    clear();
    safeDeletePtr(treePtr);
    return;
}

void STKFMM::FMMData::evaluate(std::vector<double> &srcSLValue, std::vector<double> &srcDLValue,
                               std::vector<double> &trgValue) {
    const size_t nSrc = treeDataPtr->src_coord.Dim() / 3;
    const size_t nSurf = treeDataPtr->surf_coord.Dim() / 3;
    const size_t nTrg = treeDataPtr->trg_coord.Dim() / 3;

    if (nTrg * kdimTrg != trgValue.size()) {
        printf("trg value size error for kernel %d\n", kernelChoice);
    }
    if (nSrc * kdimSL != srcSLValue.size()) {
        printf("src SL value size error for kernel %d\n", kernelChoice);
    }
    if (nSurf * kdimDL != srcDLValue.size()) {
        printf("src DL value size error for kernel %d\n", kernelChoice);
    }
    PtFMM_Evaluate(treePtr, trgValue, nTrg, &srcSLValue, &srcDLValue);
}

STKFMM::STKFMM(int multOrder_, int maxPts_, PAXIS pbc_, unsigned int kernelComb_)
    : multOrder(multOrder_), maxPts(maxPts_), pbc(pbc_), kernelComb(kernelComb_), xlow(0), xhigh(1), ylow(0), yhigh(1),
      zlow(0), zhigh(1), scaleFactor(1), xshift(0), yshift(0), zshift(0) {
    // set periodic boundary condition
    switch (pbc) {
    case PAXIS::NONE:
        periodicType = PeriodicType::NONE;
        break;
    case PAXIS::PZ:
        periodicType = PeriodicType::PZ;
        break;
    case PAXIS::PXY:
        periodicType = PeriodicType::PXY;
        break;
    case PAXIS::PXYZ:
        periodicType = PeriodicType::PXYZ;
        break;
    }
    if (pbc != PAXIS::NONE) {
        printf("to be implemented\n");
        exit(1);
    }

    poolFMM.clear();

    // parse the choice of kernels, use bitwise and
    if (kernelComb & asInteger(KERNEL::PVel)) {
        printf("enable PVel %d\n", kernelComb & asInteger(KERNEL::PVel));
        poolFMM[KERNEL::PVel] = new FMMData(KERNEL::PVel, multOrder, maxPts);
    }
    if (kernelComb & asInteger(KERNEL::PVelGrad)) {
        printf("enable PVelGrad %d\n", kernelComb & asInteger(KERNEL::PVelGrad));
        poolFMM[KERNEL::PVelGrad] = new FMMData(KERNEL::PVelGrad, multOrder, maxPts);
    }
    if (kernelComb & asInteger(KERNEL::PVelLaplacian)) {
        printf("enable PVelLaplacian %d\n", kernelComb & asInteger(KERNEL::PVelLaplacian));
        poolFMM[KERNEL::PVelLaplacian] = new FMMData(KERNEL::PVelLaplacian, multOrder, maxPts);
    }
    if (kernelComb & asInteger(KERNEL::Traction)) {
        printf("enable Traction %d\n", kernelComb & asInteger(KERNEL::Traction));
        poolFMM[KERNEL::Traction] = new FMMData(KERNEL::Traction, multOrder, maxPts);
    }
    if (kernelComb & asInteger(KERNEL::LAPPGrad)) {
        printf("enable LAPPGrad %d\n", kernelComb & asInteger(KERNEL::LAPPGrad));
        poolFMM[KERNEL::LAPPGrad] = new FMMData(KERNEL::LAPPGrad, multOrder, maxPts);
    }

#ifdef FMMDEBUG
    pvfmm::Profile::Enable(true);
#endif

    if (poolFMM.empty()) {
        printf("Error: no kernel choosed");
        exit(1);
    }

    printf("FMM Initialized\n");
}

STKFMM::~STKFMM() {
    // delete all FMMData
    for (auto &fmm : poolFMM) {
        safeDeletePtr(fmm.second);
    }
}

void STKFMM::setBox(double xlow_, double xhigh_, double ylow_, double yhigh_, double zlow_, double zhigh_) {
    xlow = xlow_;
    xhigh = xhigh_;
    ylow = ylow_;
    yhigh = yhigh_;
    zlow = zlow_;
    zhigh = zhigh_;

    // find and calculate scale & shift factor to map the box to [0,1)
    xshift = -xlow;
    yshift = -ylow;
    zshift = -zlow;
    double xlen = xhigh - xlow;
    double ylen = yhigh - ylow;
    double zlen = zhigh - zlow;
    scaleFactor = 1 / std::max(zlen, std::max(xlen, ylen));
    // new coordinate = (x+xshift)*scaleFactor, in [0,1)

    std::cout << "box x " << xlen << " box y " << ylen << " box z " << zlen << std::endl;
    std::cout << "scale factor " << scaleFactor << std::endl;

    // sanity check of box setting, ensure fitting in a cubic box [0,1)^3
    const double eps = pow(10, -12) / scaleFactor;
    switch (pbc) {
    case PAXIS::NONE:
        // for PNONE, scale max length to [0,1), all choices are valid
        break;
    case PAXIS::PZ:
        if (zlen < xlen || zlen < ylen) {
            std::cout << "periodic box size error" << std::endl;
            exit(1);
        }
        break;
    case PAXIS::PXY:
        // for PXY,PXZ,PYZ, periodic direcitons must have equal size, and larger than the third direction
        if (fabs(xlen - ylen) < eps && xlen >= zlen) {
            // correct
        } else {
            std::cout << "periodic box size error" << std::endl;
            exit(1);
        }
        break;
    case PAXIS::PXYZ:
        // for PXYZ, must be cubic
        if (fabs(xlen - ylen) < eps && fabs(xlen - zlen) < eps && fabs(ylen - zlen) < eps) {
            // correct
        } else {
            std::cout << "periodic box size error" << std::endl;
            exit(1);
        }
        break;
    }
}

void STKFMM::setupCoord(const std::vector<double> &coordIn, std::vector<double> &coord) {
    // apply scale to internal data array, without rotation
    // Set source points, with scale

    const int npts = coordIn.size() / 3;
    coord.resize(npts * 3);

    if (pbc == PAXIS::PXYZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < npts; i++) {
            coord[3 * i] = fracwrap((coordIn[3 * i] + xshift) * scaleFactor);
            coord[3 * i + 1] = fracwrap((coordIn[3 * i + 1] + yshift) * scaleFactor);
            coord[3 * i + 2] = fracwrap((coordIn[3 * i + 2] + zshift) * scaleFactor);
        }
    } else if (pbc == PAXIS::PZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < npts; i++) {
            coord[3 * i] = ((coordIn[3 * i] + xshift) * scaleFactor);
            coord[3 * i + 1] = ((coordIn[3 * i + 1] + yshift) * scaleFactor);
            coord[3 * i + 2] = fracwrap((coordIn[3 * i + 2] + zshift) * scaleFactor);
        }
    } else if (pbc == PAXIS::PXY) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < npts; i++) {
            coord[3 * i] = fracwrap((coordIn[3 * i] + xshift) * scaleFactor);
            coord[3 * i + 1] = fracwrap((coordIn[3 * i + 1] + yshift) * scaleFactor);
            coord[3 * i + 2] = ((coordIn[3 * i + 2] + zshift) * scaleFactor);
        }
    } else {
        assert(pbc == PAXIS::NONE);
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < npts; i++) {
            coord[3 * i] = ((coordIn[3 * i] + xshift) * scaleFactor);
            coord[3 * i + 1] = ((coordIn[3 * i + 1] + yshift) * scaleFactor);
            coord[3 * i + 2] = ((coordIn[3 * i + 2] + zshift) * scaleFactor);
        }
    }
    return;
}

void STKFMM::setPoints(const std::vector<double> &srcSLCoord_, const std::vector<double> &srcDLCoord_,
                       const std::vector<double> &trgCoord_) {
    int np, myrank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &myrank);

    if (!poolFMM.empty()) {
        for (auto &fmm : poolFMM) {
            printf("kernel %d \n", asInteger(fmm.second->kernelChoice));
            fmm.second->deleteTree();
        }
        printf("ALL FMM Tree Cleared\n");
    }

    // setup point coordinates
    setupCoord(srcSLCoord_, srcSLCoordInternal);
    setupCoord(srcDLCoord_, srcDLCoordInternal);
    setupCoord(trgCoord_, trgCoordInternal);
    printf("points set\n");
}

void STKFMM::setupTree(KERNEL kernel_) {
    poolFMM[kernel_]->setupTree(srcSLCoordInternal, srcDLCoordInternal, trgCoordInternal);
    printf("Coord setup for kernel %d\n", static_cast<int>(kernel_));
}

void STKFMM::evaluateFMM(std::vector<double> &srcSLValue, std::vector<double> &srcDLValue,
                         std::vector<double> &trgValue, KERNEL kernel) {

    if (poolFMM.find(kernel) == poolFMM.end()) {
        printf("Error: no such FMMData exists for kernel %d\n", static_cast<int>(kernel));
    }
    FMMData &fmm = *((*poolFMM.find(kernel)).second);

    const int nSL = srcSLCoordInternal.size() / 3;
    const int nDL = srcDLCoordInternal.size() / 3;
    const int nTrg = trgCoordInternal.size() / 3;
    srcSLValueInternal.resize(nSL * fmm.kdimSL);
    srcDLValueInternal.resize(nDL * fmm.kdimDL);
    trgValue.resize(nTrg * fmm.kdimTrg);

    // scale the source strength, SL as 1/r, DL as 1/r^2
    // SL no extra scaling
    // DL scale as scaleFactor
    srcSLValueInternal = srcSLValue;
#pragma omp parallel for
    for (int i = 0; i < nDL * fmm.kdimDL; i++) {
        srcDLValueInternal[i] = srcDLValue[i] * scaleFactor;
    }

    // run FMM
    // evaluate on internal sources with proper scaling
    trgValueInternal.resize(nTrg * fmm.kdimTrg);
    fmm.evaluate(srcSLValueInternal, srcDLValueInternal, trgValueInternal);

    // scale back according to kernel
    switch (kernel) {
    case KERNEL::PVel: {
// 1+3
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            trgValue[4 * i] = trgValueInternal[4 * i] * scaleFactor * scaleFactor; // pressure 1/r^2
            trgValue[4 * i + 1] = trgValueInternal[4 * i + 1] * scaleFactor;       // vel 1/r
            trgValue[4 * i + 2] = trgValueInternal[4 * i + 2] * scaleFactor;
            trgValue[4 * i + 3] = trgValueInternal[4 * i + 3] * scaleFactor;
        }
    } break;
    case KERNEL::PVelGrad: {
// 1+3+3+9
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            trgValue[16 * i] = trgValueInternal[16 * i] * scaleFactor * scaleFactor; // p
            for (int j = 1; j < 4; j++) {
                trgValue[16 * i + j] = trgValueInternal[16 * i + j] * scaleFactor; // vel
            }
            for (int j = 4; j < 7; j++) {
                trgValue[16 * i + j] = trgValueInternal[16 * i + j] * scaleFactor * scaleFactor * scaleFactor; // grad p
            }
            for (int j = 7; j < 16; j++) {
                trgValue[16 * i + j] = trgValueInternal[16 * i + j] * scaleFactor * scaleFactor; // grad vel
            }
        }
    } break;
    case KERNEL::Traction: {
// 9
#pragma omp parallel for
        for (int i = 0; i < 9 * nTrg; i++) {
            trgValue[i] = trgValueInternal[i] * scaleFactor * scaleFactor; // traction 1/r^2
        }
    } break;
    case KERNEL::PVelLaplacian: {
// 1+3+3
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            trgValue[7 * i] = trgValueInternal[7 * i] * scaleFactor * scaleFactor; // p
            for (int j = 1; j < 4; j++) {
                trgValue[7 * i + j] = trgValueInternal[7 * i + j] * scaleFactor; // vel
            }
            for (int j = 4; j < 7; j++) {
                trgValue[7 * i + j] =
                    trgValueInternal[7 * i + j] * scaleFactor * scaleFactor * scaleFactor; // laplacian vel
            }
        }
    } break;
    case KERNEL::LAPPGrad: {
// 1+3
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            trgValue[4 * i] = trgValueInternal[4 * i] * scaleFactor; // p, 1/r
            for (int j = 1; j < 4; j++) {
                trgValue[4 * i + j] = trgValueInternal[4 * i + j] * scaleFactor * scaleFactor; // grad p, 1/r^2
            }
        }

    } break;
    }

    return;
}

void STKFMM::showActiveKernels() {
    printf("active kernels:\n");
    if (kernelComb & asInteger(KERNEL::PVel)) {
        printf("PVel\n");
    }
    if (kernelComb & asInteger(KERNEL::PVelGrad)) {
        printf("PVelGrad\n");
    }
    if (kernelComb & asInteger(KERNEL::Traction)) {
        printf("Traction\n");
    }
    if (kernelComb & asInteger(KERNEL::PVelLaplacian)) {
        printf("PVelLaplacian\n");
    }
    if (kernelComb & asInteger(KERNEL::LAPPGrad)) {
        printf("LAPPGrad\n");
    }
}

void STKFMM::clearFMM(KERNEL kernelChoice) {
    trgValueInternal.clear();
    poolFMM[kernelChoice]->clear();
}