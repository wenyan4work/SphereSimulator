/*
 * FMMWrapper.h
 *
 *  Created on: Oct 6, 2016
 *      Author: wyan
 */

#ifndef INCLUDE_FMMWRAPPER_H_
#define INCLUDE_FMMWRAPPER_H_

// a wrapper for pvfmm
// choose kernel at compile time
#include "mpi.h"
#include "pvfmm.hpp"

#include "Eigen/Dense"

class FMM_Wrapper {
public:
    enum PAXIS {
        NONE, PXYZ, PX, PY, PZ, PXY, PXZ, PYZ
    };

    pvfmm::PtFMM matrix;
    pvfmm::PtFMM_Tree *treePtr;
    pvfmm::PtFMM_Data tree_data;

    const int mult_order;
    const int max_pts;
    const int init_depth;
    PAXIS pbc;

    double xlow, xhigh; // box
    double ylow, yhigh;
    double zlow, zhigh;
    double scaleFactor;
    double xshift, yshift, zshift;

    FMM_Wrapper(int mult_order = 8, int max_pts = 500, int init_depth = 0, PAXIS pbc_ = PAXIS::NONE);

    ~FMM_Wrapper();

    void FMM_TreeClear();

    void FMM_DataClear();

    void FMM_Evaluate(std::vector<double> &, const int, std::vector<double> *,
            std::vector<double> *surf_coordPtr = NULL);
    void FMM_UpdateTree(const std::vector<double> &, const std::vector<double> &,
            const std::vector<double> *surf_valPtr = NULL);

    void FMM_SetBox(double, double, double, double, double, double);

private:
    double *readM2LMat(const char *, const int);

    double *pm2l; // the periodizing operator

    void calcM(const pvfmm::Vector<double> &, std::vector<double> &, const std::vector<double> &);

    int pEquiv;
    int equivN;
    double scaleLEquiv;      // = 1.05;
    double pCenterLEquiv[3]; // = { -(scaleLEquiv - 1) / 2, -(scaleLEquiv - 1) / 2, -(scaleLEquiv - 1) / 2 };
    std::vector<double> M2Lsource; // the equivalent sources after the operation

    std::vector<double> pointLEquiv; // = surface(pEquiv, (double *) &(pCenterLCheck[0]), scaleLCheck, 0); // center at
                                     // 0.5,0.5,0.5, periodic box 1,1,1, scale 1.05, depth = 0
    double timeTree;
    double timeNear;
    double timeM2L;
    double timeFar;
};

template<class Real_t>
std::vector<Real_t> surface(int, Real_t *, Real_t, int);

#endif /* INCLUDE_FMMWRAPPER_H_ */
