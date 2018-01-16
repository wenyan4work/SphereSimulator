#include <cassert>
#include <cstdio>

#include "Base64.hpp"
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
    getGrid(gridPoints, gridWeights, gridValues);

    std::vector<int32_t> connect;
    std::vector<int32_t> offset;
    std::vector<uint8_t> type;
    getGridCellConnect(connect, offset, type);

    std::string contentB64; // data converted to base64 format
    contentB64.reserve(10 * order * order * (kind == KIND::LAP ? 1 : 3));
    // write a vtk unstructured grid section
    // assume filePtr in 'append' mode

    const int p = order;
    fprintf(filePtr, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 2 * p * p, 2 * p * (p + 1));

    // point data
    fprintf(filePtr, "<PointData %s=\"%s\">\n", (kind == KIND::LAP ? "Scalars" : "Vectors"), this->name.c_str());
    getBase64FromVector(gridValues, contentB64);
    fprintf(filePtr, "<DataArray Name=\"%s\">\n", this->name.c_str());
    fprintf(filePtr, "%s", contentB64);
    fprintf(filePtr, "</DataArray>\n");
    fprintf(filePtr, "</PointData>\n");
    // cell data (empty)

    // point location
    fprintf(filePtr, "<Points>\n");
    fprintf(filePtr, "<DataArray NumberOfComponents=\"3\">\n");
    contentB64.clear();
    getBase64FromVector(gridPoints, contentB64);
    fprintf(filePtr, "%s", contentB64);
    fprintf(filePtr, "</DataArray>\n");
    fprintf(filePtr, "</Points>\n");

    // cell definition
    fprintf(filePtr, "<Cells>");
    fprintf(filePtr, "<DataArray type=\"Int32\" Name=\"connectivity\">\n");
    contentB64.clear();
    getBase64FromVector(connect, contentB64);
    fprintf(filePtr, "%s\n", contentB64.c_str());
    fprintf(filePtr, "</DataArray>\n");
    fprintf(filePtr, "<DataArray type=\"Int32\" Name=\"offsets\">\n");
    contentB64.clear();
    getBase64FromVector(offset, contentB64);
    fprintf(filePtr, "%s\n", contentB64.c_str());
    fprintf(filePtr, "</DataArray>\n");
    fprintf(filePtr, "<DataArray type=\"UInt8\" Name=\"types\">\n");
    contentB64.clear();
    getBase64FromVector(type, contentB64);
    fprintf(filePtr, "%s\n", contentB64.c_str());
    fprintf(filePtr, "</DataArray>\n");

    fprintf(filePtr, "</Cells>");

    // end
    fprintf(filePtr, "</Piece>\n");
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

void SPHExp::getGridCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
                                std::vector<uint8_t> &type) const {
    const int p = order;
    // cells with north pole, 3 point for each cell
    int index = 1; // the index to the node point, starting from 1 point after the north pole
    offset.push_back(0);
    for (int k = 0; k < 2 * p + 1; k++) {
        // 3 points, 0,k,k+1
        gridCellConnect.push_back(0);
        gridCellConnect.push_back(index);
        gridCellConnect.push_back(index + 1);
        type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
        offset.push_back(offset.back() + 3);
        index++;
    }
    gridCellConnect.push_back(0);
    gridCellConnect.push_back(2 * p + 1);
    gridCellConnect.push_back(1);
    type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
    offset.push_back(offset.back() + 3);
    index++;

    // 4 cells for each cell in the center
    for (int j = 1; j < p - 1; j++) {
        for (int k = 0; k < 2 * p + 1; k++) {
            // 4 points, index, index+1, index-(2p+2), index+1 - (2p+2)
            gridCellConnect.push_back(index);
            gridCellConnect.push_back(index + 1);
            gridCellConnect.push_back(index - (2 * p + 2));
            gridCellConnect.push_back(index + 1 - (2 * p + 2));
            type.push_back(uint8_t(9)); // 9 = VTK_QUAD
            offset.push_back(offset.back() + 4);
            index++;
        }
        // last one, connect to the first one in this circle
        gridCellConnect.push_back(index);
        gridCellConnect.push_back(index - (2 * p + 1));
        gridCellConnect.push_back(index - (2 * p + 2));
        gridCellConnect.push_back(index - (2 * p + 1) - (2 * p + 2));
        type.push_back(uint8_t(9)); // 9 = VTK_QUAD
        offset.push_back(offset.back() + 4);
        index++;
    }

    // cells with south pole, 3 points for each cell
    for (int k = 0; k < 2 * p + 1; k++) {
        // 3 points, 0,k,k+1
        gridCellConnect.push_back(2 * p * p - 1);
        gridCellConnect.push_back(index);
        gridCellConnect.push_back(index + 1);
        type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
        offset.push_back(offset.back() + 3);
        index++;
    }
    gridCellConnect.push_back(2 * p * p - 1);               // south pole
    gridCellConnect.push_back(2 * p * p - 2);               // 1 before south pole
    gridCellConnect.push_back(2 * p * p - 2 - (2 * p + 1)); // 1 around the circle
    type.push_back(uint8_t(5));                             // 5= VTK_TRIANGLE

    return;
}