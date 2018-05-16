#include "SphereSPHOperator.hpp"

SphereSPHOperator::SphereSPHOperator(const std::vector<Sphere> &sphere, const std::string &name_,
                                     std::shared_ptr<STKFMM> &fmmPtr_, const double cSL_, const double cDL_,
                                     const double cTrac_)
    : spherePtr{&sphere}, 
    dimension(sphere[0].getLayer(name_).dimension), 
    name(name_), fmmPtr(fmmPtr_), cSL(cSL_),
      cDL(cDL_), cTrac(cTrac_) {
    commRcp = getMPIWORLDTCOMM();

    // tasks:
    // 1 setup the operator (fmm tree, temporary data field, etc)
    // 2 setup the right side
    // 3 setup the iterative solver
    // 4 (optional) setup the preconditioner for the operator

    // step 1 setup temporary copy and dof of sph
    setupDOF();

    // step 2 setup FMM tree and temporary space
    // FMM box should have been set out of this
    // FMM active kernels should have been set out of this
    setupFMM();
}

void SphereSPHOperator::setupDOF() {
    const auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
    sphereMapRcp = getTMAPFromLocalSize(nLocal, commRcp);

    sph.resize(nLocal);
    spectralDofIndex.resize(nLocal);
    spectralDofOffset.resize(nLocal, 0);
    gridValueDofIndex.resize(nLocal);
    gridValueDofOffset.resize(nLocal, 0);
    gridWeightDofIndex.resize(nLocal);
    gridWeightDofIndex.resize(nLocal, 0);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sph[i] = sphere[i].getLayer(name);
        spectralDofIndex[i] = sph[i].getSpectralDOF();
        gridValueDofIndex[i] = sph[i].getGridDOF();
        gridWeightDofIndex[i] = gridValueDofIndex[i] / dimension;
    }

    for (int i = 1; i < nLocal; i++) {
        spectralDofOffset[i] = spectralDofOffset[i - 1] + spectralDofIndex[i - 1];
        gridValueDofIndex[i] = gridValueDofIndex[i - 1] + gridValueDofIndex[i - 1];
        gridWeightDofIndex[i] = gridWeightDofIndex[i - 1] + gridWeightDofIndex[i - 1];
    }
    spectralDofMapRcp = getTMAPFromLocalSize(spectralDofOffset.back() + spectralDofIndex.back(), commRcp);
    gridValueDofMapRcp = getTMAPFromLocalSize(gridValueDofOffset.back() + gridValueDofIndex.back(), commRcp);
    gridWeightDofMapRcp = getTMAPFromLocalSize(gridWeightDofOffset.back() + gridWeightDofIndex.back(), commRcp);

    gridPoints.resize(3 * gridWeightDofMapRcp->getNodeNumElements());
    gridWeights.resize(gridWeightDofMapRcp->getNodeNumElements());
    gridValues.resize(gridValueDofMapRcp->getNodeNumElements());

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        std::vector<double> points;  // 3d
        std::vector<double> weights; // 1d
        std::vector<double> values;  // 1d for LAP, 3d for STK
        // returned by this contains north and south pole
        // sph[i].getGrid(points, weights, values, sphere[i].radius, sphere[i].pos);

        // remove north and south pole
        std::copy(points.cbegin() + 3, points.cend() - 3, gridPoints.begin() + 3 * gridWeightDofIndex[i]);
        std::copy(weights.cbegin() + 1, weights.cend() - 1, gridWeights.begin() + gridWeightDofIndex[i]);
        std::copy(values.cbegin() + dimension, values.cend() - dimension, gridValues.begin() + gridValueDofIndex[i]);
    }
}

void SphereSPHOperator::setupFMM() {
    assert(fmmPtr);
    if (cSL > 0 || cTrac > 0) {
        srcSLCoord = gridPoints;
    } else {
        srcSLCoord.clear();
    }

    if (cDL > 0) {
        srcDLCoord = gridPoints;
    } else {
        srcDLCoord.clear();
    }

    trgCoord = gridPoints;

    fmmPtr->setPoints(srcSLCoord, srcDLCoord, trgCoord);

    if (sph[0].kind == Shexp::KIND::LAP) {
        // laplace, SL 1d, DL 3d, trg 4d
        srcSLValue.resize(srcSLCoord.size() / 3 * 1);
        srcDLValue.resize(srcDLCoord.size() / 3 * 3);
        trgValue.resize(trgCoord.size() / 3 * 4);
        fmmPtr->setupTree(KERNEL::LAPPGrad);
    } else {
        // stokes, SL 4d, DL 9d, trg 4d (SL+DL) or 9d (Trac)
        srcSLValue.resize(srcSLCoord.size() / 3 * 4);
        srcDLValue.resize(srcDLCoord.size() / 3 * 9);
        trgValue.resize(trgCoord.size() / 3 * (cTrac == 0 ? 4 : 9));
        if (cSL > 0 || cDL > 0) {
            fmmPtr->setupTree(KERNEL::PVel);
        }
        if (cTrac > 0) {
            fmmPtr->setupTree(KERNEL::Traction);
        }
    }

    return;
}

template <class Fntr>
void SphereSPHOperator::setupRightSide(Fntr & fntr) {
    rightSideRcp = Teuchos::rcp(new TMV(spectralDofMapRcp, 1, true));

    auto rsPtr = rightSideRcp->getLocalView<Kokkos::HostSpace>();
    rightSideRcp->modify<Kokkos::HostSpace>();

    // fill entries for each sphere
    const auto & sphere = *spherePtr;
    const int nLocal = sphere.size();
    #pragma omp parallel for
    for(int i=0;i<nLocal;i++){
        std::vector<double> bvec;
        fntr(sphere[i],bvec);
    }

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSPHOperator::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, scalar_type alpha,
                              scalar_type beta) const {
    assert(mode == Teuchos::NO_TRANS);

    if (commRcp->getRank() == 0)
        printf("SphereSPHOperator Applied\n");
}