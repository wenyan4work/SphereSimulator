#include "SphereSTKSHOperator.hpp"

constexpr double pi = 3.141592653589793238462643383279;

SphereSTKSHOperator::SphereSTKSHOperator(const std::vector<Sphere> *const spherePtr, const std::string &name_,
                                         std::shared_ptr<STKFMM> &fmmPtr_, const double cIdentity_, const double cSL_,
                                         const double cTrac_, const double cLOP_)
    : spherePtr(spherePtr), name(name_), fmmPtr(fmmPtr_), cId(cIdentity_), cSL(cSL_), cTrac(cTrac_), cLOP(cLOP_) {
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
        sph.emplace_back(sphere[i].getLayer(name));
    }

    gridNumberIndex.clear();
    gridNumberIndex.resize(nLocal, 0);
    gridNumberLength.clear();
    gridNumberLength.resize(nLocal, 0);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        gridNumberLength[i] = sph[i].getGridNumber();
    }

    for (int i = 1; i < nLocal; i++) {
        gridNumberIndex[i] = gridNumberIndex[i - 1] + gridNumberLength[i - 1];
    }

    gridNumberMapRcp = getTMAPFromLocalSize(gridNumberIndex.back() + gridNumberLength.back(), commRcp);
    gridValueDofMapRcp = getTMAPFromLocalSize(gridNumberMapRcp->getNodeNumElements() * 3, commRcp);

    gridCoords.resize(3 * gridNumberMapRcp->getNodeNumElements());
    gridCoordsRelative.resize(3 * gridNumberMapRcp->getNodeNumElements());
    gridNorms.resize(3 * gridNumberMapRcp->getNodeNumElements());
    gridWeights.resize(gridNumberMapRcp->getNodeNumElements());

    TEUCHOS_ASSERT(gridCoords.size() == gridValueDofMapRcp->getNodeNumElements());
    TEUCHOS_ASSERT(gridNorms.size() == gridValueDofMapRcp->getNodeNumElements());
}

void SphereSTKSHOperator::setupFMM() {
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
            const int indexBase = gridNumberIndex[i];
            const int npts = gridNumberLength[i];
            // remove north and south pole
            std::copy(coords.cbegin() + 3, coords.cend() - 3, gridCoords.begin() + 3 * indexBase);
            std::copy(norms.cbegin() + 3, norms.cend() - 3, gridNorms.begin() + 3 * indexBase);
            std::copy(weights.cbegin() + 1, weights.cend() - 1, gridWeights.begin() + indexBase);

            // coordsrelative = norms * radius
            for (auto &v : norms) {
                v *= sph[i].radius;
            }
            std::copy(norms.cbegin() + 3, norms.cend() - 3, gridCoordsRelative.begin() + 3 * indexBase);
        }
    }

    srcSLCoord = gridCoords;
    trgCoord = gridCoords;

    const int nSL = srcSLCoord.size() / 3;
    const int nTrg = trgCoord.size() / 3;
    srcDLCoord.clear();
    srcDLValue.clear();

    fmmPtr->setPoints(nSL, srcSLCoord.data(), 0, srcDLCoord.data(), nTrg, trgCoord.data());

    // stokes, SL 4d, DL 9d, trg 4d (SL+DL) or 9d (Trac)
    srcSLValue.resize(srcSLCoord.size() / 3 * 4);
    trgValue.resize(trgCoord.size() / 3 * 9);

    // TODO: figure out a more flexible way to setup tree
    // For the moment, PVel tree and Traction tree has no double layer src points,
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
        TEUCHOS_ASSERT(length == bvec.size());
        for (int j = 0; j < length; j++) {
            rsPtr(index + j, 0) = bvec[j];
        }
    }

    return;
}

// Y := beta*Y + alpha*Op(A)*X
void SphereSTKSHOperator::apply(const TMV &X, TMV &Y, Teuchos::ETransp mode, scalar_type alpha,
                                scalar_type beta) const {
    TEUCHOS_ASSERT(mode == Teuchos::NO_TRANS);

    if (commRcp->getRank() == 0)
        printf("SphereSTKSHOperator Applied\n");

    const int nCol = X.getNumVectors();
    TEUCHOS_ASSERT(nCol == Y.getNumVectors());
    const int nRowLocal = X.getLocalLength();
    TEUCHOS_ASSERT(nRowLocal == gridValueDofMapRcp->getNodeNumElements());
    pointValues.clear();
    pointValues.resize(nRowLocal, 0);
    pointValuesApply.clear();
    pointValuesApply.resize(nRowLocal, 0);

    auto XPtr = X.getLocalView<Kokkos::HostSpace>();
    auto YPtr = Y.getLocalView<Kokkos::HostSpace>();
    Y.modify<Kokkos::HostSpace>();

    for (int c = 0; c < nCol; c++) {
#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            pointValues[i] = XPtr(i, c);
        }
        std::fill(pointValuesApply.begin(), pointValuesApply.end(), 0);

        // step 1, project out the linear space cannot be represented by spherical harmonics
        // printf("projectout\n");
        // projectNullSpace(pointValues.data());

        // step 2, run FMM
        printf("runFMM\n");
        applyP2POP(pointValues.data(), pointValuesApply.data(), cId, cSL, cTrac);
        for (auto &v : pointValuesApply) {
            printf("%lf\t", v);
        }
        printf("\n");

        // step 3, apply the rigid body operator
        printf("applyLOP\n");
        applyLOP(pointValues.data(), pointValuesApply.data(), cLOP);
        for (auto &v : pointValuesApply) {
            printf("%lf\t", v);
        }
        printf("\n");

        // step 4, store to y
        TEUCHOS_ASSERT(nRowLocal == pointValues.size());
        TEUCHOS_ASSERT(nRowLocal == pointValuesApply.size());
        TEUCHOS_ASSERT(nRowLocal == gridCoords.size());
        TEUCHOS_ASSERT(nRowLocal == gridNorms.size());
        TEUCHOS_ASSERT(nRowLocal == 3 * gridWeights.size());
#pragma omp parallel for
        for (int i = 0; i < nRowLocal; i++) {
            YPtr(i, c) = beta * YPtr(i, c) + alpha * pointValuesApply[i];
        }
    }
}

// in-place removal
void SphereSTKSHOperator::projectNullSpace(double *inPtr) const {
    const int nLocal = sph.size();
#pragma omp parallel
    {
        // temporary data space
        std::vector<double> spectralCoeff;
        std::vector<double> gridValue;
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            spectralCoeff.resize(sph[i].getSpectralDOF());
            std::fill(spectralCoeff.begin(), spectralCoeff.end(), 0);
            const int indexBase = gridNumberIndex[i];
            const int npts = gridNumberLength[i];

            gridValue.resize(3 * npts);
            std::copy(inPtr + 3 * indexBase, inPtr + 3 * indexBase + 3 * npts, gridValue.begin());

            sph[i].calcSpectralCoeff(spectralCoeff.data(), gridValue.data());
            sph[i].calcGridValue(spectralCoeff.data(), gridValue.data());

            for (int j = 0; j < npts; j++) {
                inPtr[3 * (indexBase + j) + 0] = gridValue[3 * j + 0];
                inPtr[3 * (indexBase + j) + 1] = gridValue[3 * j + 1];
                inPtr[3 * (indexBase + j) + 2] = gridValue[3 * j + 2];
            }
        }
    }
}

void SphereSTKSHOperator::applyP2POP(const double *inPtr, double *outPtr, double cIdex, double cSLex,
                                     double cTracex) const {

    std::vector<double> gridValueTemp(gridValueDofMapRcp->getNodeNumElements());
    std::copy(inPtr, inPtr + gridValueTemp.size(), gridValueTemp.data());
    cacheSHCoeff(gridValueTemp.data());

    // ex means "extra" parameter value different from stored in the matrix.
    const int nGridPts = gridNumberMapRcp->getNodeNumElements();
    const int nLocal = spherePtr->size();
    srcSLValue.clear();
    trgValue.clear();
    srcDLValue.clear();

    // SL Stokes PVel FMM
    if (fabs(cSLex) > 1e-9) {
        // step 1 setup value
        srcSLValue.resize(nGridPts * 4);
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            srcSLValue[4 * i + 0] = inPtr[3 * i + 0] * gridWeights[i] * cSLex;
            srcSLValue[4 * i + 1] = inPtr[3 * i + 1] * gridWeights[i] * cSLex;
            srcSLValue[4 * i + 2] = inPtr[3 * i + 2] * gridWeights[i] * cSLex;
            srcSLValue[4 * i + 3] = 0;
        }

        trgValue.clear();
        trgValue.resize(nGridPts * 4, 0); // trg value = PVel

        // step 2 run
        fmmPtr->clearFMM(KERNEL::PVel);
        fmmPtr->evaluateFMM(nGridPts, srcSLValue.data(), 0, srcDLValue.data(), nGridPts, trgValue.data(), KERNEL::PVel);
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            // trgValue = p,vx,vy,vz
            outPtr[3 * i + 0] += trgValue[4 * i + 1];
            outPtr[3 * i + 1] += trgValue[4 * i + 2];
            outPtr[3 * i + 2] += trgValue[4 * i + 3];
        }

        // step 3 fix trgValue with operator on self
#pragma omp parallel
        {
            std::vector<double> fmmvalue;
            std::vector<double> shvalue;
            std::vector<double> shcoeff;
            std::vector<double> gridcoordrelative;
#pragma omp for
            for (int i = 0; i < nLocal; i++) {
                const int indexBase = gridNumberIndex[i];
                const int npts = gridNumberLength[i];

                fmmvalue.resize(4 * npts);
                std::fill(fmmvalue.begin(), fmmvalue.end(), 0);

                shvalue.resize(3 * npts);
                std::fill(shvalue.begin(), shvalue.end(), 0);

                shcoeff.resize(shCoeffLength[i]);
                std::copy(shCoeffValues.cbegin() + shCoeffIndex[i],
                          shCoeffValues.cbegin() + shCoeffIndex[i] + shCoeffLength[i], shcoeff.begin());

                gridcoordrelative.resize(3 * npts);
                std::copy(gridCoordsRelative.cbegin() + 3 * indexBase,
                          gridCoordsRelative.cbegin() + 3 * indexBase + 3 * npts, gridcoordrelative.begin());

                // this is already scaled by cSLex
                fmmPtr->evaluateKernel(1, PPKERNEL::SLS2T, npts, gridcoordrelative.data(),
                                       srcSLValue.data() + 4 * indexBase, npts,
                                       gridcoordrelative.data() + 3 * indexBase, fmmvalue.data(), KERNEL::PVel);
                // this is not scaled by cSLex
                sph[i].calcSDLNF(shcoeff.data(), npts, gridcoordrelative.data(), shvalue.data(), false, true);
                for (int j = 0; j < npts; j++) {
                    outPtr[3 * (indexBase + j) + 0] += cSLex * shvalue[3 * j + 0] - fmmvalue[4 * j + 1];
                    outPtr[3 * (indexBase + j) + 1] += cSLex * shvalue[3 * j + 1] - fmmvalue[4 * j + 2];
                    outPtr[3 * (indexBase + j) + 2] += cSLex * shvalue[3 * j + 2] - fmmvalue[4 * j + 3];
                }
            }
        }
    }

    // Traction, Stokes Traction FMM
    if (fabs(cTracex) > 1e-9) {
        srcSLValue.resize(nGridPts * 4);
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            srcSLValue[4 * i] = inPtr[3 * i] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 1] = inPtr[3 * i + 1] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 2] = inPtr[3 * i + 2] * cTracex * gridWeights[i];
            srcSLValue[4 * i + 3] = 0;
        }
        trgValue.clear();
        trgValue.resize(nGridPts * 9, 0);

        fmmPtr->clearFMM(KERNEL::Traction);
        fmmPtr->evaluateFMM(nGridPts, srcSLValue.data(), 0, srcDLValue.data(), nGridPts, trgValue.data(),
                            KERNEL::Traction);
        // fix trgValue with operator on self
#pragma omp parallel
        {
            std::vector<double> fmmvalue;
            std::vector<double> shvalue;
            std::vector<double> shcoeff;
            std::vector<double> gridcoordrelative;
#pragma omp for
            for (int i = 0; i < nLocal; i++) {
                const int indexBase = gridNumberIndex[i];
                const int npts = gridNumberLength[i];

                fmmvalue.resize(9 * npts);
                std::fill(fmmvalue.begin(), fmmvalue.end(), 0);

                shvalue.resize(3 * npts);
                std::fill(shvalue.begin(), shvalue.end(), 0);

                shcoeff.resize(shCoeffLength[i]);
                std::copy(shCoeffValues.cbegin() + shCoeffIndex[i],
                          shCoeffValues.cbegin() + shCoeffIndex[i] + shCoeffLength[i], shcoeff.begin());

                gridcoordrelative.resize(3 * npts);
                std::copy(gridCoordsRelative.cbegin() + 3 * indexBase,
                          gridCoordsRelative.cbegin() + 3 * indexBase + 3 * npts, gridcoordrelative.begin());

                // this is already scaled by cTracex
                fmmPtr->evaluateKernel(1, PPKERNEL::SLS2T, npts, gridcoordrelative.data(),
                                       srcSLValue.data() + 4 * indexBase, npts,
                                       gridcoordrelative.data() + 3 * indexBase, fmmvalue.data(), KERNEL::Traction);
                // this is not scaled by cTracex
                sph[i].calcKSelf(shcoeff.data(), npts, gridcoordrelative.data(), shvalue.data(), false);

                // remove fmm self value from trgValue
                for (int j = 0; j < npts; j++) {
                    for (int k = 0; k < 9; k++) {
                        trgValue[9 * (indexBase + j) + k] -= fmmvalue[9 * j + k];
                    }
                }
                // compute velocity
                for (int j = 0; j < npts; j++) {
                    const int index = indexBase + j;
                    outPtr[3 * index + 0] += trgValue[9 * index + 0] * gridNorms[3 * index + 0] +
                                             trgValue[9 * index + 1] * gridNorms[3 * index + 1] +
                                             trgValue[9 * index + 2] * gridNorms[3 * index + 2];
                    outPtr[3 * index + 1] += trgValue[9 * index + 3] * gridNorms[3 * index + 0] +
                                             trgValue[9 * index + 4] * gridNorms[3 * index + 1] +
                                             trgValue[9 * index + 5] * gridNorms[3 * index + 2];
                    outPtr[3 * index + 2] += trgValue[9 * index + 6] * gridNorms[3 * index + 0] +
                                             trgValue[9 * index + 7] * gridNorms[3 * index + 1] +
                                             trgValue[9 * index + 8] * gridNorms[3 * index + 2];
                }

                for (int j = 0; j < npts; j++) {
                    outPtr[3 * (indexBase + j) + 0] += cTracex * shvalue[3 * j + 0];
                    outPtr[3 * (indexBase + j) + 1] += cTracex * shvalue[3 * j + 1];
                    outPtr[3 * (indexBase + j) + 2] += cTracex * shvalue[3 * j + 2];
                }
            }
        }
    }

    if (fabs(cIdex) > 1e-9) {
#pragma omp parallel for
        for (int i = 0; i < nGridPts; i++) {
            outPtr[3 * i + 0] += cIdex * inPtr[3 * i + 0];
            outPtr[3 * i + 1] += cIdex * inPtr[3 * i + 1];
            outPtr[3 * i + 2] += cIdex * inPtr[3 * i + 2];
        }
    }
}

void SphereSTKSHOperator::applyLOP(const double *inPtr, double *outPtr, double cLOPex) const {
    // apply the rigid body operator L
    auto &sphere = *spherePtr;
    const int nLocal = sphere.size();
#pragma omp parallel
    {
        // std::vector<double> temp;
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            Evec3 A(0, 0, 0);
            Evec3 B(0, 0, 0);
            const int indexBase = gridNumberIndex[i];
            const int npts = gridNumberLength[i];
            const double radius = sph[i].radius;
            for (int j = 0; j < npts; j++) {
                const int index = indexBase + j;
                Evec3 valj = Evec3(inPtr[3 * index], inPtr[3 * index + 1], inPtr[3 * index + 2]);
                A += gridWeights[index] * valj;
                B += (gridWeights[index] * radius) *
                     Evec3(gridNorms[3 * index], gridNorms[3 * index + 1], gridNorms[3 * index + 2]).cross(valj);
            }
            A *= (1 / (4 * pi * radius * radius));
            B *= (3 / (8 * pi * radius * radius * radius * radius));
            // Lj = A + B cross rj
            for (int j = 0; j < npts; j++) {
                const int index = indexBase + j;
                Evec3 valx = A + B.cross(radius * Evec3(gridNorms[3 * index], gridNorms[3 * index + 1],
                                                        gridNorms[3 * index + 2]));
                outPtr[3 * index + 0] += valx[0] * cLOPex;
                outPtr[3 * index + 1] += valx[1] * cLOPex;
                outPtr[3 * index + 2] += valx[2] * cLOPex;
            }
        }
    }
}

void SphereSTKSHOperator::cacheSHCoeff(double *gridValuesPtr) const {
    const int nLocal = sph.size();
    TEUCHOS_ASSERT(nLocal == spherePtr->size());

    shCoeffIndex.resize(nLocal, 0);
    shCoeffLength.resize(nLocal, 0);
    std::fill(shCoeffIndex.begin(), shCoeffIndex.end(), 0);
    std::fill(shCoeffLength.begin(), shCoeffLength.end(), 0);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        shCoeffLength[i] = sph[i].getSpectralDOF();
    }

    for (int i = 1; i < nLocal; i++) {
        shCoeffIndex[i] = shCoeffIndex[i - 1] + shCoeffLength[i - 1];
    }

    shCoeffValues.resize(shCoeffIndex.back() + shCoeffLength.back());
    std::fill(shCoeffValues.begin(), shCoeffValues.end(), 0);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const int indexBaseCoeff = shCoeffIndex[i];
        const int indexBaseGrid = gridNumberIndex[i];
        sph[i].calcSpectralCoeff(shCoeffValues.data() + indexBaseCoeff, gridValuesPtr + 3 * indexBaseGrid);
    }
}