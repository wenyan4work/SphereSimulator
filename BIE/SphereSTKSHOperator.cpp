#include "SphereSTKSHOperator.hpp"

constexpr double pi = 3.141592653589793238462643383279;

SphereSTKSHOperator::SphereSTKSHOperator(const std::vector<Sphere> *const spherePtr, const std::string &name_,
                                         std::shared_ptr<STKFMM> &fmmPtr_, const double cIdentity_, const double cSL_,
                                         const double cDL_, const double cTrac_, const double cLOP_)
    : spherePtr(spherePtr), name(name_), fmmPtr(fmmPtr_), cId(cIdentity_), cSL(cSL_), cDL(cDL_), cTrac(cTrac_),
      cLOP(cLOP_) {
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
    for (int i = 0; i < nLocal; i++) {
        sph.push_back(sphere[i].getLayer(name));
    }

    gridWeights.resize(nLocal);
    gridCoords.resize(nLocal * 3);
    gridNorms.resize(nLocal * 3);

    gridDofIndex.resize(nLocal, 0);
    gridDofLength.resize(nLocal);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        gridDofLength[i] = sph[i].getGridNumber();
    }

    for (int i = 1; i < nLocal; i++) {
        gridDofIndex[i] = gridDofIndex[i - 1] + gridDofLength[i - 1];
    }

    gridDofMapRcp = getTMAPFromLocalSize(gridDofIndex.back() + gridDofLength.back(), commRcp);
    gridValueDofMapRcp = getTMAPFromLocalSize(gridDofMapRcp->getNodeNumElements() * 3, commRcp);

    gridCoords.resize(3 * gridDofMapRcp->getNodeNumElements());
    gridNorms.resize(3 * gridDofMapRcp->getNodeNumElements());
    gridWeights.resize(gridDofMapRcp->getNodeNumElements());

    assert(gridCoords.size() == gridValueDofMapRcp->getNodeNumElements());
    assert(gridNorms.size() == gridValueDofMapRcp->getNodeNumElements());
}

void SphereSTKSHOperator::setupFMM() {
    assert(fmmPtr);
    const auto &sphere = *spherePtr;
    const int nLocal = sphere.size();

    // setup points
#pragma omp parallel
    {
        std::vector<double> coords;  // 3d
        std::vector<double> norms;   // 3d for STK
        std::vector<double> weights; // 1d

#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            coords.resize(sph[i].getGridNumber() * 3 + 6); // 3d
            norms.resize(sph[i].getGridNumber() * 3 + 6);  // 3d for STK
            weights.resize(sph[i].getGridNumber() + 2);    // 1d

            // returned by this contains north and south pole
            sph[i].getGridWithPole(coords, weights, sphere[i].pos, &norms);
            const int npts = gridDofLength[i];
            const int indexBase = gridDofIndex[i];
            // remove north and south pole
            std::copy(coords.cbegin() + 3, coords.cend() - 3, gridCoords.begin() + 3 * indexBase);
            std::copy(norms.cbegin() + 3, norms.cend() - 3, gridNorms.begin() + 3 * indexBase);
            std::copy(weights.cbegin() + 1, weights.cend() - 1, gridWeights.begin() + indexBase);
        }
    }

    if (cSL != 0 || cTrac != 0) {
        srcSLCoord = gridCoords;
    } else {
        srcSLCoord.clear();
    }

    if (cDL != 0) {
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
    trgValue.resize(trgCoord.size() / 3 * 9);

    fmmPtr->setupTree(KERNEL::PVel);
    fmmPtr->setupTree(KERNEL::Traction);

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
        assert(length == bvec.size());
        for (int j = 0; j < length; j++) {
            rsPtr(index + j, 0) = bvec[j];
        }
    }

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKSHOperator::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha,
                                scalar_type beta) const {
    assert(mode == Teuchos::NO_TRANS);

    if (commRcp->getRank() == 0)
        printf("SphereSTKSHOperator Applied\n");

    const int nCol = X.getNumVectors();
    assert(nCol == Y.getNumVectors());
    const int nRowLocal = X.getLocalLength();
    pointValues.clear();
    pointValues.resize(nRowLocal, 0);
    pointValuesApply.clear();
    pointValuesApply.resize(nRowLocal, 0);

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();

    for (int c = 0; c < nCol; c++) {
        for (int i = 0; i < nRowLocal; i++) {
            pointValues[i] = XPtr(i, c);
        }

        // step 1, project out the linear space cannot be represented by spherical harmonics
        printf("projectout\n");
        projectNullSpace(pointValues.data(), pointValuesApply.data());

        // step 2, run FMM
        printf("runFMM\n");
        runFMM(pointValues.data(), pointValuesApply.data(), cId, cSL, cDL, cTrac);

        // step 3, apply the rigid body operator
        printf("applyLOP\n");
        applyLOP(pointValues.data(), pointValuesApply.data());

        // step 4, store to y
        assert(nRowLocal == pointValues.size());
        assert(nRowLocal == pointValuesApply.size());
        assert(nRowLocal == gridCoords.size());
        assert(nRowLocal == gridNorms.size());
        assert(nRowLocal == 3 * gridWeights.size());
        for (int i = 0; i < nRowLocal; i++) {
            YPtr(i, c) = beta * YPtr(i, c) + alpha * pointValuesApply[i];
        }
    }
}

void SphereSTKSHOperator::projectNullSpace(const double *inPtr, double *outPtr) const {
    const int nLocal = sph.size();
#pragma omp parallel
    {
        // temporary data space
        std::vector<double> spectralCoeff;
        std::vector<double> gridValue;
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            spectralCoeff.clear();
            spectralCoeff.resize(sph[i].getSpectralDOF(), 0);
            gridValue.clear();
            gridValue.resize(sph[i].getGridDOF(), 0);
            std::copy(inPtr, inPtr + 3 * gridDofLength[i], gridValue.begin());
            sph[i].calcSpectralCoeff(spectralCoeff.data(), gridValue.data());
            sph[i].calcGridValue(spectralCoeff.data(), gridValue.data());

            const int indexBase = gridDofIndex[i];
            const int npts = gridDofLength[i];

            for (int j = 0; j < npts; j++) {
                outPtr[3 * (indexBase + j)] += gridValue[3 * j];
                outPtr[3 * (indexBase + j) + 1] += gridValue[3 * j + 1];
                outPtr[3 * (indexBase + j) + 2] += gridValue[3 * j + 2];
            }
        }
    }
}

void SphereSTKSHOperator::runFMM(const double *inPtr, double *outPtr, double cIdex, double cSLex, double cDLex,
                                 double cTracex) const {
    // ex means "extra" parameter value different from stored in the matrix.
    const int nGridPts = gridDofMapRcp->getNodeNumElements();
    const int nLocal = spherePtr->size();
    srcSLValue.resize(0);
    srcDLValue.resize(0);
    trgValue.resize(0);
    int nSL = 0, nDL = 0, nTrg = 0;

    // SL and DL, Stokes PVel FMM
    {
        // step 1 setup value
        if (cSLex != 0) {
            nSL = nGridPts;
            srcSLValue.resize(nSL * 4);
#pragma omp parallel for
            for (int i = 0; i < nSL; i++) {
                srcSLValue[4 * i] = inPtr[3 * i] * gridWeights[i] * cSLex;
                srcSLValue[4 * i + 1] = inPtr[3 * i + 1] * gridWeights[i] * cSLex;
                srcSLValue[4 * i + 2] = inPtr[3 * i + 2] * gridWeights[i] * cSLex;
                srcSLValue[4 * i + 3] = 0;
            }
        }

        if (cDLex != 0) {
            nDL = nGridPts;
            srcDLValue.resize(nDL * 9);
#pragma omp parallel for
            for (int i = 0; i < nDL; i++) {
                srcDLValue[9 * i] = gridNorms[3 * i] * inPtr[3 * i] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 1] = gridNorms[3 * i] * inPtr[3 * i + 1] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 2] = gridNorms[3 * i] * inPtr[3 * i + 2] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 3] = gridNorms[3 * i + 1] * inPtr[3 * i] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 4] = gridNorms[3 * i + 1] * inPtr[3 * i + 1] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 5] = gridNorms[3 * i + 1] * inPtr[3 * i + 2] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 6] = gridNorms[3 * i + 2] * inPtr[3 * i] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 7] = gridNorms[3 * i + 2] * inPtr[3 * i + 1] * gridWeights[i] * cDLex;
                srcDLValue[9 * i + 8] = gridNorms[3 * i + 2] * inPtr[3 * i + 2] * gridWeights[i] * cDLex;
            }
        }
        nTrg = nGridPts;
        trgValue.resize(nTrg * 3);
        // step 2 run
        if (nSL || nDL) {
            fmmPtr->evaluateFMM(nSL, srcSLValue.data(), nDL, srcDLValue.data(), nTrg, trgValue.data(), KERNEL::PVel);
#pragma omp parallel for
            for (int i = 0; i < nTrg; i++) {
                // trgValue = p,vx,vy,vz
                outPtr[3 * i] += trgValue[4 * i + 1];
                outPtr[3 * i + 1] += trgValue[4 * i + 2];
                outPtr[3 * i + 2] += trgValue[4 * i + 3];
            }
        }

        // TODO: step 3 fix trgValue with operator on self
        if (nSL) {
#pragma omp parallel for
            for (int i = 0; i < nLocal; i++) {
            }
        }

        if (nDL) {
#pragma omp parallel for
            for (int i = 0; i < nLocal; i++) {
            }
        }
    }

    // Traction, Stokes Traction FMM
    if (cTracex != 0) {
        nSL = nGridPts;
        srcSLValue.resize(nSL * 4);
#pragma omp parallel for
        for (int i = 0; i < nSL; i++) {
            srcSLValue[4 * i] = inPtr[3 * i] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 1] = inPtr[3 * i + 1] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 2] = inPtr[3 * i + 2] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 3] = 0;
        }
        nTrg = nGridPts;
        trgValue.resize(nTrg * 9);
        fmmPtr->evaluateFMM(nSL, srcSLValue.data(), 0, srcDLValue.data(), nTrg, trgValue.data(), KERNEL::Traction);
        // TODO: fix trgValue with operator on self
        // step 3 post process
#pragma omp parallel for
        for (int i = 0; i < nTrg; i++) {
            outPtr[3 * i + 0] += trgValue[9 * i + 0] * gridNorms[3 * i] + trgValue[9 * i + 1] * gridNorms[3 * i + 1] +
                                 trgValue[9 * i + 2] * gridNorms[3 * i + 2];
            outPtr[3 * i + 1] += trgValue[9 * i + 3] * gridNorms[3 * i] + trgValue[9 * i + 4] * gridNorms[3 * i + 1] +
                                 trgValue[9 * i + 5] * gridNorms[3 * i + 2];
            outPtr[3 * i + 2] += trgValue[9 * i + 6] * gridNorms[3 * i] + trgValue[9 * i + 7] * gridNorms[3 * i + 1] +
                                 trgValue[9 * i + 8] * gridNorms[3 * i + 2];
        }
    }

    if (cIdex != 0) {
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            outPtr[3 * i] += cIdex * inPtr[3 * i];
            outPtr[3 * i + 1] += cIdex * inPtr[3 * i + 1];
            outPtr[3 * i + 2] += cIdex * inPtr[3 * i + 2];
        }
    }
}
void SphereSTKSHOperator::applyLOP(const double *inPtr, double *outPtr) const {
    // apply the rigid body operator L
    auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
#pragma omp parallel
    {
        std::vector<double> temp;
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            Evec3 A(0, 0, 0);
            Evec3 B(0, 0, 0);
            const int indexBase = gridDofIndex[i];
            const int npts = gridDofLength[i];
            const double radius = sph[i].radius;
            for (int j = 0; j < npts; j++) {
                const int index = indexBase + j;
                Evec3 valj = Evec3(inPtr[3 * index], inPtr[3 * index + 1], inPtr[3 * index + 2]);
                A += gridWeights[index] * valj;
                B += gridWeights[index] * radius *
                     Evec3(gridNorms[3 * index], gridNorms[3 * index + 1], gridNorms[3 * index + 2]).cross(valj);
            }
            A *= (1 / (4 * pi * radius * radius));
            B *= (3 / (8 * pi * radius * radius * radius * radius));
            // Lj = A + B cross rj
            for (int j = 0; j < npts; j++) {
                const int index = indexBase + j;
                Evec3 valx = A + B.cross(radius * Evec3(gridNorms[3 * index], gridNorms[3 * index + 1],
                                                        gridNorms[3 * index + 2]));
                outPtr[3 * index] += valx[0] * cLOP;
                outPtr[3 * index + 1] += valx[1] * cLOP;
                outPtr[3 * index + 2] += valx[2] * cLOP;
            }
        }
    }
}