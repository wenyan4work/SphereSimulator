/*
 * FMMWrapper.cpp
 *
 *  Created on: Oct 20, 2016
 *      Author: wyan
 */

#include <cassert>

#include "FMMWrapper.h"
#include "omp.h"
#include <chrono> // for timing

PeriodicType periodicType;

// TreeClear(), optional, should be called BEFORE TreeUpdate.

// return fraction part between [0,1)
/*
 * This function is only applied in the PERIODIC DIRECTION
 *
 * The user of the library must ensure that all points are located within [0,1)
 *
 * */
inline double fracwrap(double x) {
    return x - floor(x);
}

FMM_Wrapper::FMM_Wrapper(int mult_order, int max_pts, int init_depth, PAXIS pbc_) :
        mult_order(mult_order), max_pts(max_pts), init_depth(init_depth), pbc(pbc_), xlow(0), xhigh(1), ylow(0), yhigh(
                1), zlow(0), zhigh(1), scaleFactor(1), xshift(0), yshift(0), zshift(0) {

    // set periodic boundary condition
    switch (pbc) {
        case PAXIS::NONE:
            periodicType = PeriodicType::NONE;
            break;
        case PAXIS::PXYZ:
            periodicType = PeriodicType::PXYZ;
            break;
        case PAXIS::PX:
            periodicType = PeriodicType::PZ; // use axis rotation
            break;
        case PAXIS::PY:
            periodicType = PeriodicType::PZ; // use axis rotation
            break;
        case PAXIS::PZ:
            periodicType = PeriodicType::PZ;
            break;
        case PAXIS::PXY:
            periodicType = PeriodicType::PXY;
            break;
        case PAXIS::PXZ:
            periodicType = PeriodicType::PXY; // use axis rotation
            break;
        case PAXIS::PYZ:
            periodicType = PeriodicType::PXY; // use axis rotation
            break;
    }
    pm2l = nullptr;
    if (pbc != NONE) {
        if (mult_order != (mult_order / 2) * 2 || mult_order < 6 || mult_order > 16) {
            printf("periodic M2L data available only for p=6,8,10,12,14,16\n");
        } else if (pbc == PAXIS::PXYZ) {
            switch (mult_order) {
                case 6:
                    pm2l = readM2LMat("M2LStokes3D3Dp6", 6);
                    break;
                case 8:
                    pm2l = readM2LMat("M2LStokes3D3Dp8", 8);
                    break;
                case 10:
                    pm2l = readM2LMat("M2LStokes3D3Dp10", 10);
                    break;
                case 12:
                    pm2l = readM2LMat("M2LStokes3D3Dp12", 12);
                    break;
                case 14:
                    pm2l = readM2LMat("M2LStokes3D3Dp14", 14);
                    break;
                case 16:
                    pm2l = readM2LMat("M2LStokes3D3Dp16", 16);
                    break;
                default:
                    std::cout << "no m2l data at corresponding p, exit now" << std::endl;
                    exit(1);
                    break;
            }
        } else if (pbc == PAXIS::PX || pbc == PAXIS::PY || pbc == PAXIS::PZ) {
            switch (mult_order) {
                case 6:
                    pm2l = readM2LMat("M2LStokes1D3Dp6", 6);
                    break;
                case 8:
                    pm2l = readM2LMat("M2LStokes1D3Dp8", 8);
                    break;
                case 10:
                    pm2l = readM2LMat("M2LStokes1D3Dp10", 10);
                    break;
                case 12:
                    pm2l = readM2LMat("M2LStokes1D3Dp12", 12);
                    break;
                case 14:
                    pm2l = readM2LMat("M2LStokes1D3Dp14", 14);
                    break;
                case 16:
                    pm2l = readM2LMat("M2LStokes1D3Dp16", 16);
                    break;
                default:
                    std::cout << "no m2l data at corresponding p, exit now" << std::endl;
                    exit(1);
                    break;
            }
        } else if (pbc == PAXIS::PXY || pbc == PAXIS::PXZ || pbc == PAXIS::PYZ) {
            switch (mult_order) {
                case 6:
                    pm2l = readM2LMat("M2LStokes2D3Dp6", 6);
                    break;
                case 8:
                    pm2l = readM2LMat("M2LStokes2D3Dp8", 8);
                    break;
                case 10:
                    pm2l = readM2LMat("M2LStokes2D3Dp10", 10);
                    break;
                case 12:
                    pm2l = readM2LMat("M2LStokes2D3Dp12", 12);
                    break;
                case 14:
                    pm2l = readM2LMat("M2LStokes2D3Dp14", 14);
                    break;
                case 16:
                    pm2l = readM2LMat("M2LStokes2D3Dp16", 16);
                    break;
                default:
                    std::cout << "no m2l data at corresponding p, exit now" << std::endl;
                    exit(1);
                    break;
            }
        }

        this->pEquiv = mult_order; // (8-1)^2*6 + 2 points

        this->scaleLEquiv = RAD1; // RAD1 = 2.95 defined in pvfmm_common.h
        this->pCenterLEquiv[0] = -(scaleLEquiv - 1) / 2;
        this->pCenterLEquiv[1] = -(scaleLEquiv - 1) / 2;
        this->pCenterLEquiv[2] = -(scaleLEquiv - 1) / 2;

        pointLEquiv = surface(pEquiv, (double *) &(pCenterLEquiv[0]), scaleLEquiv, 0); // center at 0.5,0.5,0.5, periodic box 1,1,1, scale 1.05, depth = 0

        equivN = 6 * (pEquiv - 1) * (pEquiv - 1) + 2;
    }
    const pvfmm::Kernel<double> &kernel_fn = pvfmm::StokesKernel<double>::velocity();

    MPI_Comm comm = MPI_COMM_WORLD;
    // treePtr = new pvfmm::PtFMM_Tree(comm);
    treePtr = nullptr;
    matrix.Initialize(mult_order, comm, &kernel_fn);
#ifdef FMMDEBUG
    pvfmm::Profile::Enable(true);
#endif

    printf("FMM Initialized\n");
}

FMM_Wrapper::~FMM_Wrapper() {
    FMM_TreeClear();
    if (pm2l != nullptr) {
        delete[] pm2l;
        pm2l = nullptr;
    }
}

void FMM_Wrapper::FMM_TreeClear() {
    FMM_DataClear();

    if (treePtr != nullptr) {
        treePtr->ClearFMMData();
        delete treePtr;
        treePtr = nullptr;
    }
}

void FMM_Wrapper::FMM_DataClear() {
    if (treePtr != nullptr) {
        treePtr->ClearFMMData();
    }
}

void FMM_Wrapper::FMM_SetBox(double xlow_, double xhigh_, double ylow_, double yhigh_, double zlow_, double zhigh_) {
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
// new coordinate = (x+xshift)*scaleFactor, in (0,1)

    std::cout << "box x" << xlen << "box y" << ylen << "box z" << zlen << std::endl;
    std::cout << "scale factor" << scaleFactor << std::endl;

// validate box setting, ensure fitting in a cubic box [0,1]^3
    const double eps = pow(10, -10) / scaleFactor;
    switch (pbc) {
        case PAXIS::NONE:
            // for PNONE, scale max length to (0,1), all choices are valid
            break;
        case PAXIS::PX:
            // for PX,PY,PZ, max must be the periodic direction
            if (xlen < ylen || xlen < zlen) {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
            break;
        case PAXIS::PY:
            if (ylen < xlen || ylen < zlen) {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
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
            } else {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
            break;
        case PAXIS::PXZ:
            if (fabs(xlen - zlen) < eps && xlen >= ylen) {
            } else {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
            break;
        case PAXIS::PYZ:
            if (fabs(zlen - ylen) < eps && zlen >= xlen) {
            } else {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
            break;
        case PAXIS::PXYZ:
            // for PXYZ, must be cubic
            if (fabs(xlen - ylen) < eps && fabs(xlen - ylen) < eps && fabs(xlen - zlen) < eps) {
            } else {
                std::cout << "periodic box size error" << std::endl;
                exit(1);
            }
            break;
    }
}

void FMM_Wrapper::FMM_UpdateTree(const std::vector<double> &src_coord, const std::vector<double> &trg_coord,
        const std::vector<double> *surf_coordPtr) {
    int np, myrank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &myrank);

    FMM_TreeClear();

    treePtr = new pvfmm::PtFMM_Tree(comm);

    tree_data.dim = 3;
    tree_data.max_depth = 15;
    tree_data.max_pts = max_pts;

// Set source points, with scale
//	tree_data.src_coord = src_coord;
    const int nsrc = src_coord.size() / 3;
    tree_data.src_coord.Resize(nsrc * 3);
    if (pbc == PAXIS::PY) {
// rotate y axis to z axis to use the z 1d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i + 1] = ((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = fracwrap((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i] = ((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PXZ) {
// rotate y axis to z axis to use the xy 2d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i + 1] = fracwrap((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = ((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i] = fracwrap((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PX) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i + 2] = fracwrap((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i] = ((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = ((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PYZ) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i + 2] = ((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i] = fracwrap((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = fracwrap((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }

    } else if (pbc == PAXIS::PXYZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i] = fracwrap((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = fracwrap((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = fracwrap((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i] = ((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = ((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = fracwrap((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PXY) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i] = fracwrap((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = fracwrap((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = ((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else {
        assert(pbc == PAXIS::NONE);
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            tree_data.src_coord[3 * i] = ((src_coord[3 * i] + xshift) * scaleFactor);
            tree_data.src_coord[3 * i + 1] = ((src_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.src_coord[3 * i + 2] = ((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.src_coord[3 * i] << tree_data.src_coord[3 * i + 1] << tree_data.src_coord[3 * i + 2]
            << std::endl;
#endif
        }
    }
    if (surf_coordPtr != nullptr) {
        // set to NULL. currently no support for surf source
        const int nsurf = 0;
        tree_data.surf_coord.Resize(nsurf * 3);
        surf_coordPtr = nullptr;
    }

// Set target points.
// use the same rotation and periodic wrap as source

    const int ntrg = trg_coord.size() / 3;
    tree_data.trg_coord.Resize(ntrg * 3);
    if (pbc == PAXIS::PY) {
// rotate y axis to z axis to use the z 1d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i + 1] = ((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = fracwrap((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i] = ((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PXZ) {
// rotate y axis to z axis to use the xy 2d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i + 1] = fracwrap((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = ((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i] = fracwrap((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PX) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i + 2] = fracwrap((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i] = ((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = ((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PYZ) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i + 2] = ((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i] = fracwrap((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = fracwrap((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }

    } else if (pbc == PAXIS::PXYZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i] = fracwrap((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = fracwrap((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = fracwrap((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PZ) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i] = ((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = ((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = fracwrap((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else if (pbc == PAXIS::PXY) {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i] = fracwrap((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = fracwrap((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = ((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    } else {
        assert(pbc == PAXIS::NONE);
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < ntrg; i++) {
            tree_data.trg_coord[3 * i] = ((trg_coord[3 * i] + xshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 1] = ((trg_coord[3 * i + 1] + yshift) * scaleFactor);
            tree_data.trg_coord[3 * i + 2] = ((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
            std::cout << tree_data.trg_coord[3 * i] << tree_data.trg_coord[3 * i + 1] << tree_data.trg_coord[3 * i + 2]
            << std::endl;
#endif
        }
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    tree_data.pt_coord = tree_data.trg_coord;

    treePtr->Initialize(&tree_data);
    bool adap = true;

    treePtr->InitFMM_Tree(adap, pbc == NONE ? pvfmm::FreeSpace : pvfmm::Periodic);
    treePtr->SetupFMM(&matrix);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "Tree construction time:" << duration / 1e6 << std::endl;
    timeTree = duration / 1e6;
#ifdef FMMDEBUG
    std::cout << "SetupFMM Complete" << std::endl;
#endif
}

void FMM_Wrapper::FMM_Evaluate(std::vector<double> &trg_val, const int n_trg, std::vector<double> *src_val,
        std::vector<double> *surf_valPtr) {
// in place rotate of src_val;
    if (src_val == nullptr) {
        printf("Error, no source value\n");
        return;
    }
    if (surf_valPtr != nullptr) {
        printf("Error, no srcval not fully implemented\n");
        return;
    }
    const int nsrc = src_val->size() / 3;
    if (pbc == PAXIS::PY || pbc == PAXIS::PXZ) {
// rotate y axis to z axis to use the z 1d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            double x = (*src_val)[3 * i];
            double y = (*src_val)[3 * i + 1];
            double z = (*src_val)[3 * i + 2];
            (*src_val)[3 * i] = z;
            (*src_val)[3 * i + 1] = x;
            (*src_val)[3 * i + 2] = y;
        }
    } else if (pbc == PAXIS::PX || pbc == PAXIS::PYZ) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            double x = (*src_val)[3 * i];
            double y = (*src_val)[3 * i + 1];
            double z = (*src_val)[3 * i + 2];
            (*src_val)[3 * i] = y;
            (*src_val)[3 * i + 1] = z;
            (*src_val)[3 * i + 2] = x;
        }
    } else {
        // no rotate
    }
#ifdef FMMTIMING
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    PtFMM_Evaluate(treePtr, trg_val, n_trg, src_val, surf_valPtr);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << "PtFMM Evaluation time:" << duration / 1e6 << std::endl;
    timeNear = duration / 1e6;
#else
    PtFMM_Evaluate(treePtr, trg_val, n_trg, src_val, surf_valPtr);
#endif

#ifdef FMMDEBUG
    std::cout << "before pxyz" << trg_val[0] << std::endl;
    std::cout << trg_val[1] << std::endl;
    std::cout << trg_val[2] << std::endl;
#endif
    if (pbc != NONE) {
#ifdef FMMTIMING
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        calcM(tree_data.trg_coord, trg_val, *src_val);
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "Periodic M2L Evaluation time:" << duration / 1e6 << std::endl;
#else
        calcM(tree_data.trg_coord, trg_val, *src_val);
#endif
    }

// scale and rotate back
    if (pbc == PAXIS::PY || pbc == PAXIS::PXZ) {
// rotate y axis to z axis to use the z 1d periodic data
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            double x = (*src_val)[3 * i];
            double y = (*src_val)[3 * i + 1];
            double z = (*src_val)[3 * i + 2];
            (*src_val)[3 * i] = y;
            (*src_val)[3 * i + 1] = z;
            (*src_val)[3 * i + 2] = x;
        }
#pragma omp parallel for
        for (size_t i = 0; i < n_trg; i++) {
            double x = trg_val[3 * i];
            double y = trg_val[3 * i + 1];
            double z = trg_val[3 * i + 2];
            trg_val[3 * i] = y * scaleFactor;
            trg_val[3 * i + 1] = z * scaleFactor;
            trg_val[3 * i + 2] = x * scaleFactor;
        }
    } else if (pbc == PAXIS::PX || pbc == PAXIS::PYZ) {
// rotate x axis to z axis
#pragma omp parallel for
        for (size_t i = 0; i < nsrc; i++) {
            double x = (*src_val)[3 * i];
            double y = (*src_val)[3 * i + 1];
            double z = (*src_val)[3 * i + 2];
            (*src_val)[3 * i] = z;
            (*src_val)[3 * i + 1] = x;
            (*src_val)[3 * i + 2] = y;
        }
#pragma omp parallel for
        for (size_t i = 0; i < n_trg; i++) {
            double x = trg_val[3 * i];
            double y = trg_val[3 * i + 1];
            double z = trg_val[3 * i + 2];
            trg_val[3 * i] = z * scaleFactor;
            trg_val[3 * i + 1] = x * scaleFactor;
            trg_val[3 * i + 2] = y * scaleFactor;
        }
    } else {
// no rotate
#pragma omp parallel for
        for (size_t i = 0; i < n_trg; i++) {
            trg_val[3 * i] *= scaleFactor;
            trg_val[3 * i + 1] *= scaleFactor;
            trg_val[3 * i + 2] *= scaleFactor;
        }
    }
}

double *FMM_Wrapper::readM2LMat(const char *fname, const int p) {
    const int size = 3 * (6 * (p - 1) * (p - 1) + 2);
    double *fdata = new double[size * size];

    char *pvfmm_dir = getenv("PVFMM_DIR");
    std::stringstream st;
    st << pvfmm_dir;
    st << "/pdata/";
    st << fname;
    FILE *fin = fopen(st.str().c_str(), "r");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int iread, jread;
            double fread;
            fscanf(fin, "%d %d %lf\n", &iread, &jread, &fread);
            if (i != iread || j != jread) {
                printf("read ij error \n");
            }
            fdata[i * size + j] = fread;
        }
    }

    fclose(fin);
    return fdata;
}

/**
 * \brief Returns the coordinates of points on the surface of a cube.
 * \param[in] p Number of points on an edge of the cube is (n+1)
 * \param[in] c Coordinates to the centre of the cube (3D array).
 * \param[in] alpha Scaling factor for the size of the cube.
 * \param[in] depth Depth of the cube in the octree.
 * \return Vector with coordinates of points on the surface of the cube in the
 * format [x0 y0 z0 x1 y1 z1 .... ].
 */

template<class Real_t>
std::vector<Real_t> surface(int p, Real_t *c, Real_t alpha, int depth) {
    size_t n_ = (6 * (p - 1) * (p - 1) + 2); // Total number of points.

    std::vector<Real_t> coord(n_ * 3);
    coord[0] = coord[1] = coord[2] = -1.0;
    size_t cnt = 1;
    for (int i = 0; i < p - 1; i++)
        for (int j = 0; j < p - 1; j++) {
            coord[cnt * 3] = -1.0;
            coord[cnt * 3 + 1] = (2.0 * (i + 1) - p + 1) / (p - 1);
            coord[cnt * 3 + 2] = (2.0 * j - p + 1) / (p - 1);
            cnt++;
        }
    for (int i = 0; i < p - 1; i++)
        for (int j = 0; j < p - 1; j++) {
            coord[cnt * 3] = (2.0 * i - p + 1) / (p - 1);
            coord[cnt * 3 + 1] = -1.0;
            coord[cnt * 3 + 2] = (2.0 * (j + 1) - p + 1) / (p - 1);
            cnt++;
        }
    for (int i = 0; i < p - 1; i++)
        for (int j = 0; j < p - 1; j++) {
            coord[cnt * 3] = (2.0 * (i + 1) - p + 1) / (p - 1);
            coord[cnt * 3 + 1] = (2.0 * j - p + 1) / (p - 1);
            coord[cnt * 3 + 2] = -1.0;
            cnt++;
        }
    for (size_t i = 0; i < (n_ / 2) * 3; i++)
        coord[cnt * 3 + i] = -coord[i];

    Real_t r = 0.5 * pow(0.5, depth);
    Real_t b = alpha * r;
    for (size_t i = 0; i < n_; i++) {
        coord[i * 3 + 0] = (coord[i * 3 + 0] + 1.0) * b + c[0];
        coord[i * 3 + 1] = (coord[i * 3 + 1] + 1.0) * b + c[1];
        coord[i * 3 + 2] = (coord[i * 3 + 2] + 1.0) * b + c[2];
    }
    return coord;
}

// fast inverse sqrt for double
inline double invsqrt(double number) {
    double y = number;
    double x2 = y * 0.5;
    std::int64_t i = *(std::int64_t *) &y;
// The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double *) &i;
    y = y * (1.5 - (x2 * y * y)); // 1st iteration
    y = y * (1.5 - (x2 * y * y)); //
    y = y * (1.5 - (x2 * y * y)); //
    y = y * (1.5 - (x2 * y * y)); //
// iteration 2-4 are necessary for p>6
    return y;
}

inline double invsqrtSimple(double number) {
    return 1 / sqrt(number);
}

void FMM_Wrapper::calcM(const pvfmm::Vector<double> &trg_coord, std::vector<double> &trg_value,
        const std::vector<double> &src_value) {
    pvfmm::Vector<double> v = treePtr->RootNode()->FMMData()->upward_equiv; // the value calculated by pvfmm
// make a copy to do correction

    assert(v.Dim() == 3 * this->equivN);
// add to trg_value
    const int n_trg = trg_coord.Dim() / 3;
    const double pi = 3.1415926535897932384626433;
    const int n_src = src_value.size() / 3;

    int M = 3 * equivN;
    int N = 3 * equivN; // checkN = equivN in this code.
    M2Lsource.resize(3 * equivN);
    assert(M2Lsource.size() == v.Dim());

    {
#ifdef FMMTIMING
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif
#pragma omp parallel for
        for (int i = 0; i < M; i++) {
            double temp = 0;
            //#pragma unroll 4
            for (int j = 0; j < N; j++) {
                temp += pm2l[i * N + j] * v[j];
            }
            M2Lsource[i] = temp;
        }
#ifdef FMMTIMING
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "Apply M2L operator time:" << duration / 1e6 << std::endl;
        timeM2L = duration / 1e6;
#endif
    }
    const double factor8pi = 1 / (8 * (double) 3.1415926535897932384626433);

    {
#ifdef FMMTIMING
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
#endif

#pragma omp parallel for
        for (int i = 0; i < n_trg; i++) {
            const double tx = trg_coord[3 * i];
            const double ty = trg_coord[3 * i + 1];
            const double tz = trg_coord[3 * i + 2];
            double trgValueX = 0;
            double trgValueY = 0;
            double trgValueZ = 0;
// rely on compiler for auto-vectorization
#pragma omp simd
            for (int p = 0; p < equivN; p++) {
                const double lx = pointLEquiv[3 * p];
                const double ly = pointLEquiv[3 * p + 1];
                const double lz = pointLEquiv[3 * p + 2];
                const double fx = M2Lsource[3 * p];
                const double fy = M2Lsource[3 * p + 1];
                const double fz = M2Lsource[3 * p + 2];
                const double rx = (tx - lx);
                const double ry = (ty - ly);
                const double rz = (tz - lz);
                const double rnorm2 = rx * rx + ry * ry + rz * rz;
                const double rinv = invsqrt(rnorm2);
                const double rinv3 = rinv * rinv * rinv;
                const double commonFac = (rx * fx + ry * fy + rz * fz);
                trgValueX += fx * rinv + commonFac * rx * rinv3;
                trgValueY += fy * rinv + commonFac * ry * rinv3;
                trgValueZ += fz * rinv + commonFac * rz * rinv3;
            }
            // do not forget the 8 pi factor, applied in the calculation step
            // the scaleFactor is applied later
            trg_value[3 * i] += trgValueX * factor8pi;
            trg_value[3 * i + 1] += trgValueY * factor8pi;
            trg_value[3 * i + 2] += trgValueZ * factor8pi;
        }
#ifdef FMMTIMING
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        std::cout << "L to point Evaluation time:" << duration / 1e6 << std::endl;
        timeFar = duration / 1e6;
#endif
    }

#ifdef FMMTIMING
    std::cout << "cost percentage = " << (timeFar) * 100 / (timeTree + timeNear) << std::endl;
    timeFar = 0;
    timeTree = 0;
    timeNear = 0;
    timeM2L = 0;
#endif

}
