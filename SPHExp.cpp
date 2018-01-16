#include <cassert>
#include <cstdio>

#include "Gauss_Legendre_Nodes_and_Weights.hpp"
#include "SPHExp.hpp"

constexpr double pi = 3.1415926535897932384626433;

// constructor
SPHExp::SPHExp(const KIND kind_, const std::string &name_, const int order_, const Equatn orientation_)
    : kind(kind_), order(order_), dimension((kind_ == KIND::LAP) ? 1 : 3), name(name_),
      orientation(orientation_) { // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL

    spectralCoeff.resize(dimension * (order + 1) * (order + 1)); // coefficients. d*(p+1)^2 elements
}

// copy constructor
SPHExp::SPHExp(const SPHExp &other)
    : kind(other.kind), order(other.order), dimension((other.kind == KIND::LAP) ? 1 : 3), name(other.name),
      orientation(other.orientation) {
    spectralCoeff = other.spectralCoeff;
}

SPHExp::SPHExp(SPHExp &&other)
    : kind(other.kind), order(other.order), dimension((other.kind == KIND::LAP) ? 1 : 3), name(other.name),
      orientation(other.orientation) {
    spectralCoeff = std::move(other.spectralCoeff);
}

void SPHExp::dumpSpectralValues(const std::string &filename) const {
    std::string kindName;
    switch (kind) {
    case KIND::LAP:
        kindName = "Lap";
        break;
    case KIND::STK:
        kindName = "Stk";
        break;
    }

    if (filename.empty()) {
        // default, output to screen
        printf("kind %s, p %8d, name %s\n", kindName.c_str(), order, name.c_str());
        for (int n = 0; n <= order; n++) {
            for (int m = -n; m <= n; m++) {
                if (dimension == 1) {
                    printf("n %3d, m %3d, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)]);
                } else if (dimension == 3) {
                    printf("n %3d, m %3d, %8f, %8f, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)],
                           spectralCoeff[COEFFINDEX(1, n, m)], spectralCoeff[COEFFINDEX(2, n, m)]);
                }
            }
        }
    } else {
        // dump to file if filename is set
        FILE *fdump = fopen(filename.c_str(), "w");
        fprintf(fdump, "kind %s, p %8d, name %s\n", kindName.c_str(), order, name.c_str());
        for (int n = 0; n <= order; n++) {
            for (int m = -n; m <= n; m++) {
                if (dimension == 1) {
                    fprintf(fdump, "n %3d, m %3d, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)]);
                } else if (dimension == 3) {
                    fprintf(fdump, "n %3d, m %3d, %8f, %8f, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)],
                            spectralCoeff[COEFFINDEX(1, n, m)], spectralCoeff[COEFFINDEX(2, n, m)]);
                }
            }
        }
        fclose(fdump);
    }
}

void SPHExp::dumpVTK(FILE *const filePtr) const {
    assert(filePtr != nullptr);

    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    std::vector<double> gridValues;

    getGrid(gridPoints,gridWeights,gridValues);

    // write a vtk unstructured grid section
    // assume filePtr in 'append' mode
}

void SPHExp::getGrid(std::vector<double> &gridPoints, std::vector<double> &gridWeights,
                     std::vector<double> &gridValues) const {
    /*
        point order: 0 at northpole, then 2p+2 points per circle. the last at south pole
    */
    const int p = order;
    gridPoints.resize(2 * p * p * 3);
    gridWeights.resize(2 * p * p * ((kind == KIND::LAP) ? 1 : 3));
    gridValues.resize(2 * p * p);

    std::vector<double> nodesGL; // cos thetaj = tj
    std::vector<double> weightsGL;

    Gauss_Legendre_Nodes_and_Weights(p + 1, nodesGL, weightsGL);

    // calculate grid coordinate without rotation
    // in total 1 + (2p+2)*(p-1) +1 = 2p^2 points
    // north pole
    int index = 0;
    gridPoints[3 * index] = 0;
    gridPoints[3 * index + 1] = 0;
    gridPoints[3 * index + 2] = 1;
    gridWeights[index] = weightsGL[0];
    index++;
    // between north and south pole
    for (int j = 1; j < p; j++) {
#pragma omp simd
        for (int k = 0; k <= 2 * p + 1; k++) {
            // sin thetaj cos phik, sin thetaj sin phik, cos thetak
            const double costhetaj = nodesGL[j];
            const double phik = 2 * pi * k / (2 * p + 2);
            double sinthetaj = sqrt(1 - costhetaj * costhetaj);
            gridPoints[3 * index] = sinthetaj * cos(phik);
            gridPoints[3 * index + 1] = sinthetaj * sin(phik);
            gridPoints[3 * index + 2] = costhetaj;
            gridWeights[index] = weightsGL[j];
            index++;
        }
    }
    // south pole
    gridPoints[3 * index] = 0;
    gridPoints[3 * index + 1] = 0;
    gridPoints[3 * index + 2] = -1;
    gridWeights[index] = weightsGL[p];
    assert(index == 2 * p * p - 1);

    // fill value with zero
    // TODO: also rotate the values if vector-valued (Stokes)
    std::fill(gridValues.begin(), gridValues.end(), 0);

    // rotation with quaternion for each point
    const EAmat3 rotation = orientation.toRotationMatrix(); // an aligned temporary rotation matrix for all points
    for (int i = 0; i < 2 * p * p; i++) {
        EAvec3 pointVec(gridPoints[3 * i], gridPoints[3 * i + 1], gridPoints[3 * i + 2]); // aligned temporary object
        Eigen::Map<Eigen::Vector3d>(gridPoints.data() + 3 * i) = rotation * pointVec;
    }

    return;
}