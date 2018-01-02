#ifndef FMMINTERFACE_H
#define FMMINTERFACE_H

#include "pvfmm.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>

void stfmm3dparttarg_pvfmm(int *ier, int *iprec, int *nsrc, double *src_coord, int *ifsingle, double *srcSL,
                           int *ifdouble, double *srcDL, double *srcDV, int *ifpot, double *srcPot, double *srcPre,
                           int *ifgrad, double *srcPotGrad, int *ntrg, double *trg_coord, int *ifpottrg, double *trgPot,
                           double *trgPre, int *ifgradtrg, double *trgPotGrad);

// subroutine stfmm3dpartself (ier,iprec,nparts,source, ifsingle,sigma sl,ifdouble,sigma dl,sigma dv,
// ifpot,pot,pre,ifgrad,grad)
void stfmm3dpartself_pvfmm(int *ier, int *iprec, int *nsrc, double *src_coord, int *ifsingle, double *srcSL,
                           int *ifdouble, double *srcDL, double *srcDV, int *ifpot, double *srcPot, double *srcPre,
                           int *ifgrad, double *srcPotGrad);

#endif