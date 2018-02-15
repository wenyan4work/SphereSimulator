#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "SPHExp.hpp"
#include "Util/Base64.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"
#include "Util/sctl.hpp"

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

int SPHExp::dumpVTK(std::ofstream &file, const std::array<double, 3> &coordBase) const {
    // indexBase is the index of the first grid point

    // this must be called in single thread
    // change to using cpp to avoid copy from string to c_str()

    assert(file.is_open());
    assert(file.good());

    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    std::vector<double> gridValues;
    getGrid(gridPoints, gridWeights, gridValues);
    const int nPts = gridWeights.size();
    assert(gridPoints.size() == nPts * 3);

    for (int i = 0; i < nPts; i++) {
        gridPoints[3 * i] += coordBase[0];
        gridPoints[3 * i + 1] += coordBase[1];
        gridPoints[3 * i + 2] += coordBase[2];
    }

    // for debug
    printf("%d,%d,%d\n", gridPoints.size(), gridWeights.size(), gridValues.size());
    for (const auto &v : gridPoints) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    for (const auto &v : gridWeights) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::string contentB64; // data converted to base64 format
    contentB64.reserve(10 * order * order * (kind == KIND::LAP ? 1 : 3));
    contentB64.clear();

    std::vector<char> temp{'a', 'b', 'c', 'd', 'e', 'f'};
    B64Converter::getBase64FromVector(temp, contentB64);
    for (const auto &v : contentB64) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::vector<int32_t> connect;
    std::vector<int32_t> offset;
    std::vector<uint8_t> type;
    getGridCellConnect(connect, offset, type);

    // for (auto &v : connect) {
    //     v += indexBase;
    // }

    // for debug
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

    // write a vtk unstructured grid section
    // assume file in 'append' mode

    const int p = order;
    file << "<Piece NumberOfPoints=\"" << 2 * p * p + 4 * p + 4 << "\" NumberOfCells=\"" << (2 * p + 2) * (p + 2)
         << "\">\n";

    // point data
    file << "<PointData Scalars=\"scalars\">\n";
    contentB64.clear();
    B64Converter::getBase64FromVector(gridValues, contentB64);
    file << "<DataArray Name=\"" << this->name << "\" type=\"Float64\" NumberOfComponents=\""
         << ((kind == KIND::LAP) ? 1 : 3) << "\" format=\"binary\">\n";
    file << contentB64 << "\n";
    file << "</DataArray>\n";

    contentB64.clear();
    B64Converter::getBase64FromVector(gridWeights, contentB64);
    file << "<DataArray Name=\""
         << "weights"
         << "\" type=\"Float64\" "
         << " format=\"binary\">\n";
    file << contentB64 << "\n";
    file << "</DataArray>\n";
    file << "</PointData>\n";

    // cell data (empty)

    // point location
    file << "<Points>\n";
    file << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"binary\">\n";
    contentB64.clear();
    B64Converter::getBase64FromVector(gridPoints, contentB64);
    file << contentB64 << "\n";
    file << "</DataArray>\n";
    file << "</Points>\n";

    // cell definition
    file << "<Cells>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
    contentB64.clear();
    B64Converter::getBase64FromVector(connect, contentB64);
    file << contentB64 << "\n";
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
    contentB64.clear();
    B64Converter::getBase64FromVector(offset, contentB64);
    file << contentB64 << "\n";
    file << "</DataArray>\n";
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n";
    contentB64.clear();
    B64Converter::getBase64FromVector(type, contentB64);
    file << contentB64 << "\n";
    file << "</DataArray>\n";
    file << "</Cells>\n";

    // end
    file << "</Piece>" << std::endl; // flush

    return nPts;
}

void SPHExp::getGrid(std::vector<double> &gridPoints, std::vector<double> &gridWeights,
                     std::vector<double> &gridValues) const {
    /*
        point order: 0 at northpole, then 2p+2 points per circle. the last at south pole
        the north and south pole are not included in the nodesGL of Gauss-Legendre nodes.
        add those two points with weight = 0 manually.
        total point = (p+1)(2p+2) + north/south pole = 2p^2+4p+4
    */
    const int p = order;
    gridPoints.resize((2 * p * p + 4 * p + 4) * 3);
    gridValues.resize((2 * p * p + 4 * p + 4) * (kind == KIND::LAP ? 1 : 3));
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

    // fill value with zero
    // TODO: convert spectral values to grid values
    // TODO: also rotate the values if vector-valued (Stokes)
    std::fill(gridValues.begin(), gridValues.end(), 0);

    // rotation with quaternion for each point
    const EAmat3 rotation = orientation.toRotationMatrix(); // an aligned temporary rotation matrix for all points
    for (int i = 0; i < 2 * p * p + 2; i++) {
        EAvec3 pointVec(gridPoints[3 * i], gridPoints[3 * i + 1], gridPoints[3 * i + 2]); // aligned temporary object
        Eigen::Map<Eigen::Vector3d>(gridPoints.data() + 3 * i) = rotation * pointVec;
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

void SPHExp::calcGridValues(std::vector<double> &val, const std::vector<double> &coeff) const {
    typedef double Real;
    long Ncoeff = (order + 1) * (order + 2);
    long Ngrid = (order + 1) * (2 * order + 2);
    long dof = coeff.size() / Ncoeff;
    assert(coeff.size() == dof * Ncoeff);
    if (val.size() != dof * Ngrid)
        val.resize(dof * Ngrid);
    const sctl::Vector<Real> coeff_(
        coeff.size(), sctl::Ptr2Itr<Real>(coeff.size() ? (Real *)coeff.data() : nullptr, coeff.size()), false);
    sctl::Vector<Real> val_(val.size(), sctl::Ptr2Itr<Real>(val.size() ? val.data() : nullptr, val.size()), false);
    sctl::SphericalHarmonics<Real>::SHC2Grid(coeff_, sctl::SHCArrange::ROW_MAJOR, order, order + 1, 2 * order + 2,
                                             &val_);
}

void SPHExp::calcSpectralValues(std::vector<double> &coeff, const std::vector<double> &val) const {
    typedef double Real;
    long Ncoeff = (order + 1) * (order + 2);
    long Ngrid = (order + 1) * (2 * order + 2);
    long dof = val.size() / Ngrid;
    assert(val.size() == dof * Ngrid);
    if (coeff.size() != dof * Ncoeff)
        coeff.resize(dof * Ncoeff);
    sctl::Vector<Real> coeff_(coeff.size(), sctl::Ptr2Itr<Real>(coeff.size() ? coeff.data() : nullptr, coeff.size()),
                              false);
    const sctl::Vector<Real> val_(val.size(),
                                  sctl::Ptr2Itr<Real>(val.size() ? (Real *)val.data() : nullptr, val.size()), false);
    sctl::SphericalHarmonics<Real>::Grid2SHC(val_, order + 1, 2 * order + 2, order, coeff_,
                                             sctl::SHCArrange::ROW_MAJOR);
}
