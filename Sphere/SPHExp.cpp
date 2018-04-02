#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "SPHExp.hpp"
#include "Util/Base64.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"
#include "Util/IOHelper.hpp"
#include "Util/sctl.hpp"

constexpr double pi = 3.1415926535897932384626433;

// constructor
SPHExp::SPHExp(const KIND kind_, const std::string &name_, const int order_, const Equatn orientation_)
    : kind(kind_), order(order_), name(name_),
      orientation(orientation_) { // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL
    gridValues.resize(getGridDOF() +
                      getDimension() * 2); // add extra two points of north and south pole in the beginning and the end
}

// move constructor
SPHExp::SPHExp(SPHExp &&other)
    : kind(other.kind), order(other.order), name(other.name), orientation(other.orientation) {
    gridValues = std::move(other.gridValues);
}

SPHExp &SPHExp::operator=(SPHExp &&other) {
    kind = other.kind;
    order = other.order;
    name = other.name;
    orientation = other.orientation;
    gridValues = std::move(other.gridValues);
}

void SPHExp::dumpSpectralValues(const double *const spectralCoeff,
                                const std::string &filename = std::string("")) const {
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
                if (getDimension() == 1) {
                    printf("n %3d, m %3d, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)]);
                } else if (getDimension() == 3) {
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
                if (getDimension() == 1) {
                    fprintf(fdump, "n %3d, m %3d, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)]);
                } else if (getDimension() == 3) {
                    fprintf(fdump, "n %3d, m %3d, %8f, %8f, %8f\n", n, m, spectralCoeff[COEFFINDEX(0, n, m)],
                            spectralCoeff[COEFFINDEX(1, n, m)], spectralCoeff[COEFFINDEX(2, n, m)]);
                }
            }
        }
        fclose(fdump);
    }
}

// data and weight always in Float64
int SPHExp::writeVTU(std::ofstream &file, const double &radius, const Evec3 &coordBase) const {
    // indexBase is the index of the first grid point

    // this must be called in single thread
    // change to using cpp to avoid copy from string to c_str()

    assert(file.is_open());
    assert(file.good());

    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    getGrid(gridPoints, gridWeights, radius, coordBase);
    const int nPts = gridWeights.size();
    assert(gridPoints.size() == nPts * 3);

// for debug
#ifdef DEBUGVTU
    printf("%lu,%lu,%lu\n", gridPoints.size(), gridWeights.size(), gridValues.size());
    for (const auto &v : gridPoints) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    for (const auto &v : gridWeights) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
#endif

    std::string contentB64; // data converted to base64 format
    contentB64.reserve(10 * order * order * (kind == KIND::LAP ? 1 : 3));
    contentB64.clear();

    std::vector<int32_t> connect;
    std::vector<int32_t> offset;
    std::vector<uint8_t> type;
    getGridCellConnect(connect, offset, type);

// for debug
#ifdef DEBUGVTU
    for (const auto &v : connect) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    for (const auto &v : offset) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    for (const auto &v : type) {
        int temp = v;
        std::cout << temp << " ";
    }
    std::cout << std::endl;
#endif

    // write a vtk unstructured grid section
    // assume file in 'append' mode

    const int p = order;
    file << "<Piece NumberOfPoints=\"" << 2 * p * p + 4 * p + 4 << "\" NumberOfCells=\"" << (2 * p + 2) * (p + 2)
         << "\">\n";

    // point data
    file << "<PointData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(gridValues, this->name, ((kind == KIND::LAP) ? 1 : 3), file);
    IOHelper::writeDataArrayBase64(gridWeights, "weights", 1, file);
    file << "</PointData>\n";

    // cell data (empty)

    // point location
    file << "<Points>\n";
    IOHelper::writeDataArrayBase64(gridPoints, "points", 3, file);
    file << "</Points>\n";

    // cell definition
    file << "<Cells>\n";
    IOHelper::writeDataArrayBase64(connect, "connectivity", 1, file);
    IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
    IOHelper::writeDataArrayBase64(type, "types", 1, file);
    file << "</Cells>\n";

    // end
    file << "</Piece>" << std::endl; // flush

    return nPts;
}

void SPHExp::getGrid(std::vector<double> &gridPoints, std::vector<double> &gridWeights, const double &radius,
                     const Evec3 &coordBase) const {
    /*
        point order: 0 at northpole, then 2p+2 points per circle. the last at south pole
        the north and south pole are not included in the nodesGL of Gauss-Legendre nodes.
        add those two points with weight = 0 manually.
        total point = (p+1)(2p+2) + north/south pole = 2p^2+4p+4
    */
    const int p = order;
    gridPoints.resize((2 * p * p + 4 * p + 4) * 3);
    gridWeights.resize(2 * p * p + 4 * p + 4);

    std::vector<double> nodesGL; // cos thetaj = tj
    std::vector<double> weightsGL;

    Gauss_Legendre_Nodes_and_Weights(p + 1, nodesGL, weightsGL); // p+1 points, excluding the two poles

    // calculate grid coordinate without rotation
    // north pole
    int index = 0;
    gridPoints[3 * index] = 0;
    gridPoints[3 * index + 1] = 0;
    gridPoints[3 * index + 2] = 1;
    gridWeights[index] = 0;
    index++;
    // between north and south pole
    // from north pole (1) to south pole (-1), picking the points from nodesGL in reversed order
    for (int j = 0; j < p + 1; j++) {
#pragma omp simd
        for (int k = 0; k <= 2 * p + 1; k++) {
            // sin thetaj cos phik, sin thetaj sin phik, cos thetak
            const double costhetaj = nodesGL[p - j];
            const double phik = 2 * pi * k / (2 * p + 2);
            double sinthetaj = sqrt(1 - costhetaj * costhetaj);
            gridPoints[3 * index] = sinthetaj * cos(phik);
            gridPoints[3 * index + 1] = sinthetaj * sin(phik);
            gridPoints[3 * index + 2] = costhetaj;
            gridWeights[index] = weightsGL[p - j];
            index++;
        }
    }
    // south pole
    gridPoints[3 * index] = 0;
    gridPoints[3 * index + 1] = 0;
    gridPoints[3 * index + 2] = -1;
    gridWeights[index] = 0;
    index++;
    assert(index == 2 * p * p + 4 * p + 4);

    // rotation with quaternion for each point coordinate
    for (int i = 0; i < 2 * p * p + 4 * p + 4; i++) {
        EAvec3 pointVec(gridPoints[3 * i], gridPoints[3 * i + 1], gridPoints[3 * i + 2]); // aligned temporary object
        Eigen::Map<Evec3>(gridPoints.data() + 3 * i) = radius * (orientation * pointVec) + coordBase;
    }

    // rotation with quaternion for each point vectorial value
    if (kind == KIND::STK) {
        for (int i = 0; i < 2 * p * p + 4 * p + 4; i++) {
            EAvec3 pointVec(gridValues[3 * i], gridValues[3 * i + 1],
                            gridValues[3 * i + 2]); // aligned temporary object
        }
    }

    return;
}

void SPHExp::getGridCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
                                std::vector<uint8_t> &type) const {
    // offset gives the END position of each cell

    const int p = order;
    // for p=0 two points on equator. Define 4 3-node cells

    // cells with north pole, 3 point for each cell
    int index = 0; // the index to the node point, starting from 1 point after the north pole
    for (int k = 0; k < 2 * p + 1; k++) {
        // 3 points, 0,k,k+1
        index++;
        gridCellConnect.push_back(0);
        gridCellConnect.push_back(index);
        gridCellConnect.push_back(index + 1);
        type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
        offset.push_back(gridCellConnect.size());
    }
    index++;
    gridCellConnect.push_back(0);
    gridCellConnect.push_back(2 * p + 2);
    gridCellConnect.push_back(1);
    type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
    offset.push_back(gridCellConnect.size());

    // 4 cells for each cell in the center
    for (int j = 1; j < p + 1; j++) {
        for (int k = 0; k < 2 * p + 1; k++) {
            // 4 points, index, index+1, index-(2p+2), index+1 - (2p+2)
            index++;
            gridCellConnect.push_back(index);
            gridCellConnect.push_back(index + 1);
            gridCellConnect.push_back(index + 1 - (2 * p + 2));
            gridCellConnect.push_back(index - (2 * p + 2));
            type.push_back(uint8_t(9)); // 9 = VTK_QUAD
            offset.push_back(gridCellConnect.size());
        }
        // last one, connect to the first one in this circle
        index++;
        gridCellConnect.push_back(index);
        gridCellConnect.push_back(index - (2 * p + 1));
        gridCellConnect.push_back(index - (2 * p + 1) - (2 * p + 2));
        gridCellConnect.push_back(index - (2 * p + 2));
        type.push_back(uint8_t(9)); // 9 = VTK_QUAD
        offset.push_back(gridCellConnect.size());
    }

    // cells with south pole, 3 points for each cell
    gridCellConnect.push_back(2 * p * p + 4 * p + 3);               // index for south pole
    gridCellConnect.push_back(2 * p * p + 4 * p + 2);               // 1 before south pole
    gridCellConnect.push_back(2 * p * p + 4 * p + 2 - (2 * p + 1)); // 1 around the circle
    type.push_back(uint8_t(5));                                     // 5= VTK_TRIANGLE
    offset.push_back(gridCellConnect.size());
    for (int k = 0; k < 2 * p + 1; k++) {
        // 3 points, k,k+1, southpole
        index--;
        gridCellConnect.push_back(2 * p * p + 4 * p + 3); // index for south pole
        gridCellConnect.push_back(index + 1);
        gridCellConnect.push_back(index);
        type.push_back(uint8_t(5)); // 5= VTK_TRIANGLE
        offset.push_back(gridCellConnect.size());
    }

    return;
}

// not including the poles
void SPHExp::calcGridValue(double *valPtr, double *coeffPtr) const {
    typedef double Real;
    // ROW_MAJOR = (order+1)(order+2)/2 complex numbers, need (order+1)(order+2) doubles per component
    const int Ncoeff = (order + 1) * (order + 2);
    // grid values are all real numbers, need (order+1)*(2order+2) doubles, excluding the north and south pole
    const int Ngrid = (order + 1) * (2 * order + 2);

    const int dimension = getDimension();

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> val(Ngrid * dimension, sctl::Ptr2Itr<Real>(valPtr ? valPtr : nullptr, Ngrid * dimension), false);

    if (kind == KIND::LAP) {
        sctl::SphericalHarmonics<Real>::SHC2Grid(coeff, sctl::SHCArrange::ROW_MAJOR, order, order + 1, 2 * order + 2,
                                                 &val);
    } else {
        sctl::SphericalHarmonics<Real>::VecSHC2Grid(coeff, sctl::SHCArrange::ROW_MAJOR, order, order + 1, 2 * order + 2,
                                                    val);
    }
}
void SPHExp::calcGridValuePole(double *valPtr, double *coeffPtr) const {

    typedef double Real;
    const int Ncoeff = (order + 1) * (order + 2);
    const int Ngrid = 2; // north and south pole
    const int dimension = getDimension();

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> val(Ngrid * dimension, sctl::Ptr2Itr<Real>(valPtr ? valPtr : nullptr, Ngrid * dimension), false);

    sctl::Vector<Real> cos_theta_phi(4);
    cos_theta_phi[0] = 1;  // north pole, cos theta
    cos_theta_phi[1] = 0;  // north pole, phi
    cos_theta_phi[2] = -1; // south pole, cos theta
    cos_theta_phi[3] = 0;  // south pole, phi

    if (kind == KIND::LAP) {
        sctl::SphericalHarmonics<Real>::SHCEval(coeff, sctl::SHCArrange::ROW_MAJOR, order, cos_theta_phi, val);
    } else {
        sctl::SphericalHarmonics<Real>::VecSHCEval(coeff, sctl::SHCArrange::ROW_MAJOR, order, cos_theta_phi, val);
    }
}

void SPHExp::calcSpectralValue(double *coeffPtr, double *valPtr) const {
    typedef double Real;
    const int Ncoeff = (order + 1) * (order + 2);
    const int Ngrid = (order + 1) * (2 * order + 2);
    const int dimension = getDimension();

    sctl::Vector<Real> coeff(Ncoeff * dimension, sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension),
                             false);
    const sctl::Vector<Real> val(Ngrid * dimension, sctl::Ptr2Itr<Real>(valPtr ? valPtr : nullptr, Ngrid * dimension),
                                 false);

    if (kind == KIND::LAP) {
        sctl::SphericalHarmonics<Real>::Grid2SHC(val, order + 1, 2 * order + 2, order, coeff,
                                                 sctl::SHCArrange::ROW_MAJOR);
    } else {
        sctl::SphericalHarmonics<Real>::Grid2VecSHC(val, order + 1, 2 * order + 2, order, coeff,
                                                    sctl::SHCArrange::ROW_MAJOR);
    }
}
