#include "SphereSTKSHOperator.hpp"

SphereSTKSHOperator::SphereSTKSHOperator(const std::vector<Sphere> &sphere, const std::string &name_,
                                         std::shared_ptr<STKFMM> &fmmPtr_, const double cIdentity_, const double cSL_,
                                         const double cDL_, const double cTrac_, const double cLOP_)
    : spherePtr{&sphere}, dimension(sphere[0].getLayer(name_).getGridDOF()), name(name_), fmmPtr(fmmPtr_), cSL(cSL_),
      cDL(cDL_), cTrac(cTrac_), cLOP(cLOP_) {
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

void SphereSTKSHOperator::setupDOF() {
    const auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
    sphereMapRcp = getTMAPFromLocalSize(nLocal, commRcp);

    sph.resize(nLocal);
    gridValueDofLength.resize(nLocal);
    gridWeights.resize(nLocal);
    gridCoords.resize(nLocal * 3);
    gridNorms.resize(nLocal * 3);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sph[i] = sphere[i].getLayer(name);
        gridValueDofLength[i] = sph[i].getGridDOF();
        gridDofLength[i] = sph[i].getGridNumber();
    }

    gridValueDofIndex.resize(nLocal, 0);
    gridDofIndex.resize(nLocal, 0);

    for (int i = 1; i < nLocal; i++) {
        gridValueDofIndex[i] = gridValueDofIndex[i - 1] + gridValueDofLength[i - 1];
        gridDofIndex[i] = gridDofIndex[i - 1] + gridDofLength[i - 1];
    }

    gridValueDofMapRcp = getTMAPFromLocalSize(gridValueDofIndex.back() + gridValueDofLength.back(), commRcp);
    gridDofMapRcp = getTMAPFromLocalSize(gridDofIndex.back() + gridDofLength.back(), commRcp);

    gridValues.resize(gridValueDofMapRcp->getNodeNumElements());
    gridCoords.resize(3 * gridDofMapRcp->getNodeNumElements());
    gridNorms.resize(3 * gridDofMapRcp->getNodeNumElements());
    gridWeights.resize(3 * gridDofMapRcp->getNodeNumElements());

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        std::vector<double> coords(sph[i].getGridNumber() * 3 + 6); // 3d
        std::vector<double> weights(sph[i].getGridNumber() + 2);    // 1d
        std::vector<double> norms(sph[i].getGridNumber() * 3 + 6);  // 3d for STK
        // returned by this contains north and south pole
        // sph[i].getGrid(points, weights, values, sphere[i].radius, sphere[i].pos);
        sph[i].getGridWithPole(coords, weights, (*spherePtr)[i].pos, &norms);

        // remove north and south pole
        std::copy(coords.cbegin() + 3, coords.cend() - 3, gridCoords.begin() + 3 * gridDofIndex[i]);
        std::copy(norms.cbegin() + 3, norms.cend() - 3, gridNorms.begin() + 3 * gridDofIndex[i]);
        std::copy(weights.cbegin() + 1, weights.cend() - 1, gridWeights.begin() + gridDofIndex[i]);
    }
}

void SphereSTKSHOperator::setupFMM() {
    assert(fmmPtr);
    const auto &sphere = *spherePtr;
    const int nLocal = sphere.size();

    // setup points
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        std::vector<double> coords(sph[i].getGridNumber() * 3 + 6); // 3d
        std::vector<double> norms(sph[i].getGridNumber() * 3 + 6);  // 3d for STK
        std::vector<double> weights(sph[i].getGridNumber() + 2);    // 1d

        // returned by this contains north and south pole
        // sph[i].getGrid(points, weights, values, sphere[i].radius, sphere[i].pos);
        sph[i].getGridWithPole(coords, weights, sphere[i].pos, &norms);
        const int npts = sph[i].getGridNumber();
        // remove north and south pole
        std::copy(coords.cbegin() + 3, coords.cend() - 3, gridCoords.begin() + 3 * npts);
        std::copy(norms.cbegin() + 3, norms.cend() - 3, gridNorms.begin() + 3 * npts);
        std::copy(weights.cbegin() + 1, weights.cend() - 1, gridWeights.begin() + npts);
    }

    if (cSL > 0 || cTrac > 0) {
        srcSLCoord = gridCoords;
    } else {
        srcSLCoord.clear();
    }

    if (cDL > 0) {
        srcDLCoord = gridCoords;
    } else {
        srcDLCoord.clear();
    }

    trgCoord = gridCoords;

    const int nSL = srcSLCoord.size() / 3;
    const int nDL = srcDLCoord.size() / 3;
    const int nTrg = trgCoord.size() / 3;

    fmmPtr->setPoints(nSL, srcSLCoord.data(), nDL, srcDLCoord.data(), nTrg, trgCoord.data());

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

    return;
}

template <class Fntr>
void SphereSTKSHOperator::setupRightSide(Fntr &fntr) {
    // right side has the same number of dofs for each sphere
    rightSideRcp = Teuchos::rcp(new TMV(spectralDofMapRcp, 1, true));

    auto rsPtr = rightSideRcp->getLocalView<Kokkos::HostSpace>();
    rightSideRcp->modify<Kokkos::HostSpace>();

    // fill entries for each sphere
    const auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        std::vector<double> bvec;
        fntr(sphere[i], bvec);

        const int index = gridValueDofIndex[i];
        const int length = gridValueDofLength[i];
        assert(length = bvec.size());
        for (int j = 0; j < length; j++) {
            rsPtr(index + j, 0) = bvec[j];
        }
    }

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKSHOperator::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, scalar_type alpha,
                                scalar_type beta) const {
    assert(mode == Teuchos::NO_TRANS);

    if (commRcp->getRank() == 0)
        printf("SphereSTKSHOperator Applied\n");

    const int nCol = X.getNumVectors();
    assert(nCol == Y.getNumVectors());
    const int nRowLocal = X.getLocalLength();
    pointValues.resize(nRowLocal);

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();

    for (int c = 0; c < nCol; c++) {
        for (int i = 0; i < nRowLocal; i++) {
            pointValues[i] = XPtr(i, c);
        }

        // step 1, project out the linear space cannot be represented by spherical harmonics
        projectNullSpace();

        // step 2, run FMM
        runFMM();

        // step 3, apply the rigid body operator
    }
}

void SphereSTKSHOperator::projectNullSpace() const {
    const int nLocal = sph.size();
#pragma omp parallel
    {

        std::vector<double> spectralCoeff;
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            spectralCoeff.resize(sph[i].getSpectralDOF(), 0);
            sph[i].calcSpectralCoeff(spectralCoeff.data(), pointValues.data() + gridValueDofIndex[i]);
            sph[i].calcGridValue(spectralCoeff.data(), pointValues.data() + gridValueDofIndex[i]);
        }
    }
}

void SphereSTKSHOperator::runFMM() const {
    const int nGridPts = gridDofMapRcp->getNodeNumElements();
    srcSLValue.resize(0);
    srcDLValue.resize(0);
    trgValue.resize(0);
    int nSL = 0, nDL = 0, nTrg = 0;

    // step 1 setup value
    if (cTrac == 0) {
        if (cSL != 0) {
            nSL = nGridPts;
            srcSLValue.resize(nGridPts * 4);
#pragma omp parallel for
            for (int i = 0; i < nGridPts; i++) {
                srcSLValue[4 * i] = pointValues[3 * i] * cSL * gridWeights[i];
                srcSLValue[4 * i + 1] = pointValues[3 * i + 1] * cSL * gridWeights[i];
                srcSLValue[4 * i + 2] = pointValues[3 * i + 2] * cSL * gridWeights[i];
                srcSLValue[4 * i + 3] = 0;
            }
        }

        if (cDL != 0) {
            nDL = nGridPts;
            srcDLValue.resize(nGridPts * 9);
#pragma omp parallel for
            for (int i = 0; i < nGridPts; i++) {
                srcDLValue[9 * i] = cDL * gridNorms[3 * i] * pointValues[3 * i] * gridWeights[i];
                srcDLValue[9 * i + 1] = cDL * gridNorms[3 * i] * pointValues[3 * i + 1] * gridWeights[i];
                srcDLValue[9 * i + 2] = cDL * gridNorms[3 * i] * pointValues[3 * i + 2] * gridWeights[i];
                srcDLValue[9 * i + 3] = cDL * gridNorms[3 * i + 1] * pointValues[3 * i] * gridWeights[i];
                srcDLValue[9 * i + 4] = cDL * gridNorms[3 * i + 1] * pointValues[3 * i + 1] * gridWeights[i];
                srcDLValue[9 * i + 5] = cDL * gridNorms[3 * i + 1] * pointValues[3 * i + 2] * gridWeights[i];
                srcDLValue[9 * i + 6] = cDL * gridNorms[3 * i + 2] * pointValues[3 * i] * gridWeights[i];
                srcDLValue[9 * i + 7] = cDL * gridNorms[3 * i + 2] * pointValues[3 * i + 1] * gridWeights[i];
                srcDLValue[9 * i + 8] = cDL * gridNorms[3 * i + 2] * pointValues[3 * i + 2] * gridWeights[i];
            }
        }
        nTrg = nGridPts;
        trgValue.resize(nGridPts * 3);
    } else {
        nSL = nGridPts;
        srcSLValue.resize(nGridPts * 4);
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            srcSLValue[4 * i] = pointValues[3 * i] * cSL * gridWeights[i];
            srcSLValue[4 * i + 1] = pointValues[3 * i + 1] * cSL * gridWeights[i];
            srcSLValue[4 * i + 2] = pointValues[3 * i + 2] * cSL * gridWeights[i];
            srcSLValue[4 * i + 3] = 0;
        }
        nTrg = nGridPts;
        trgValue.resize(nGridPts * 9);
    }

    // step 2 run
    if (cTrac == 0) {
        fmmPtr->evaluateFMM(nSL, srcSLValue.data(), nDL, srcDLValue.data(), nTrg, trgValue.data(), KERNEL::PVel);
    } else {
        fmmPtr->evaluateFMM(nSL, srcSLValue.data(), nDL, srcDLValue.data(), nTrg, trgValue.data(), KERNEL::Traction);
    }

    // step 3 post process
    if (cTrac == 0) {
        pointValues = trgValue;
    } else {
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            pointValues[3 * i + 0] = trgValue[9 * i + 0] * gridNorms[3 * i] +
                                     trgValue[9 * i + 1] * gridNorms[3 * i + 1] +
                                     trgValue[9 * i + 2] * gridNorms[3 * i + 2];
            pointValues[3 * i + 1] = trgValue[9 * i + 3] * gridNorms[3 * i] +
                                     trgValue[9 * i + 4] * gridNorms[3 * i + 1] +
                                     trgValue[9 * i + 5] * gridNorms[3 * i + 2];
            pointValues[3 * i + 2] = trgValue[9 * i + 6] * gridNorms[3 * i] +
                                     trgValue[9 * i + 7] * gridNorms[3 * i + 1] +
                                     trgValue[9 * i + 8] * gridNorms[3 * i + 2];
        }
    }
}