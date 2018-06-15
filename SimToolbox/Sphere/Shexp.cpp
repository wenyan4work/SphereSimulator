#include "Shexp.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>

#include "Util/Base64.hpp"
#include "Util/Gauss_Legendre_Nodes_and_Weights.hpp"
#include "Util/IOHelper.hpp"
#include "sctl/sctl.hpp"

constexpr double pi = 3.1415926535897932384626433;

// constructor
Shexp::Shexp(const KIND kind_, const std::string &name_, const int order_, const double &radius_,
             const Equatn orientation_)
    : kind(kind_), order(order_), name(name_), radius(radius_),
      orientation(orientation_) { // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL
    gridValue.resize(getGridDOF());
}

// move constructor
Shexp::Shexp(Shexp &&other) : kind(other.kind), order(other.order), name(other.name), orientation(other.orientation) {
    gridValue = std::move(other.gridValue);
}

Shexp &Shexp::operator=(Shexp &&other) {
    kind = other.kind;
    order = other.order;
    name = other.name;
    orientation = other.orientation;
    gridValue = std::move(other.gridValue);
}

void Shexp::dumpGridValue(const std::string &filename) const {
    std::string kindName;
    switch (kind) {
    case KIND::LAP:
        kindName = "Lap";
        break;
    case KIND::STK:
        kindName = "Stk";
        break;
    }

    const int N = order + 1;     // theta 0 to pi
    const int M = 2 * order + 2; // phi 0 to 2pi
    if (filename.empty()) { // default, output to screen printf("kind %s, p %8d, name %s\n", kindName.c_str(), order,
                            // name.c_str());
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                if (kind == KIND::LAP) {
                    printf("n %3d, m %3d, %8f\n", n, m, gridValue[m + n * M]);
                } else if (kind == KIND::STK) {
                    const int index = m + n * M;
                    printf("n %3d, m %3d, %8f, %8f, %8f\n", n, m, gridValue[3 * index], gridValue[3 * index + 1],
                           gridValue[3 * index + 2]);
                }
            }
        }
    } else {
        // dump to file if filename is set
        FILE *fdump = fopen(filename.c_str(), "w");
        fprintf(fdump, "kind %s, p %8d, name %s\n", kindName.c_str(), order, name.c_str());
        for (int n = 0; n < N; n++) {
            for (int m = 0; m < M; m++) {
                if (kind == KIND::LAP) {
                    fprintf(fdump, "n %3d, m %3d, %8f\n", n, m, gridValue[m + n * M]);
                } else if (kind == KIND::STK) {
                    const int index = m + n * M;
                    fprintf(fdump, "n %3d, m %3d, %8f, %8f, %8f\n", n, m, gridValue[3 * index],
                            gridValue[3 * index + 1], gridValue[3 * index + 2]);
                }
            }
        }
        fclose(fdump);
    }
}

// data and weight always in Float64
int Shexp::writeVTU(std::ofstream &file, const Evec3 &coordBase) const {
    // indexBase is the index of the first grid point

    // this must be called in single thread
    // change to using cpp to avoid copy from string to c_str()

    assert(file.is_open());
    assert(file.good());

    std::vector<double> gridPoints;
    std::vector<double> gridWeights;
    getGridWithPole(gridPoints, gridWeights, coordBase);
    const int nPts = gridWeights.size();
    assert(gridPoints.size() == nPts * 3);

// for debug
#ifdef DEBUGVTU
    printf("%lu,%lu,%lu\n", gridPoints.size(), gridWeights.size(), gridValue.size());
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
    getGridWithPoleCellConnect(connect, offset, type);

    std::vector<double> gridValueWithPole(getGridDOF() + 2 * getDimension());
    std::vector<double> spectralValues(getSpectralDOF());

    calcSpectralCoeff(spectralValues.data());
    double poleValues[6] = {0, 0, 0, 0, 0, 0};
    calcPoleValue(spectralValues.data(), poleValues);
    // put the pole values in the beginning and end of the array
    if (kind == KIND::LAP) {
        gridValueWithPole[0] = poleValues[0];
        std::copy(gridValue.cbegin(), gridValue.cend(), gridValueWithPole.begin() + 1);
        gridValueWithPole.back() = poleValues[1];
    } else {
        gridValueWithPole[0] = poleValues[0];
        gridValueWithPole[1] = poleValues[1];
        gridValueWithPole[2] = poleValues[2];
        std::copy(gridValue.cbegin(), gridValue.cend(), gridValueWithPole.begin() + 3);
        gridValueWithPole[3 + gridValue.size()] = poleValues[3];
        gridValueWithPole[4 + gridValue.size()] = poleValues[4];
        gridValueWithPole[5 + gridValue.size()] = poleValues[5];
        // printf("north pole: %lf,%lf,%lf\n", poleValues[0], poleValues[1], poleValues[2]);
        // printf("south pole: %lf,%lf,%lf\n", poleValues[3], poleValues[4], poleValues[5]);
    }

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
    IOHelper::writeDataArrayBase64(gridValueWithPole, this->name, ((kind == KIND::LAP) ? 1 : 3), file);
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

void Shexp::getGridWithPole(std::vector<double> &gridPoints, std::vector<double> &gridWeights, const Evec3 &coordBase,
                            std::vector<double> *gridNormPtr) const {
    /*
        point order: 0 at northpole, then 2p+2 points per circle. the last at south pole
        the north and south pole are not included in the nodesGL of Gauss-Legendre nodes.
        add those two points with weight = 0 manually.
        total point = (p+1)(2p+2) + north/south pole = 2p^2+4p+4
    */
    const int Ngrid = getGridNumber() + 2; // with poles
    gridPoints.resize(3 * Ngrid);
    gridWeights.resize(Ngrid);

    std::vector<double> nodesGL; // cos thetaj = tj
    std::vector<double> weightsGL;

    Gauss_Legendre_Nodes_and_Weights(order + 1, nodesGL, weightsGL); // p+1 points, excluding the two poles

    // calculate grid coordinate without rotation
    // north pole
    gridPoints[0] = 0;
    gridPoints[1] = 0;
    gridPoints[2] = 1;
    gridWeights[0] = 0;

    // between north and south pole
    // from north pole (1) to south pole (-1), picking the points from nodesGL in reversed order
    const int p = order;
    const double weightfactor = radius * radius * 2 * pi / (2 * p + 2);
    for (int j = 0; j < p + 1; j++) {
        for (int k = 0; k < 2 * p + 2; k++) {
            // sin thetaj cos phik, sin thetaj sin phik, cos thetak
            const double costhetaj = nodesGL[p - j];
            const double phik = 2 * pi * k / (2 * p + 2);
            const double sinthetaj = sqrt(1 - costhetaj * costhetaj);
            const int index = (j * (2 * p + 2)) + k + 1;
            gridPoints[3 * index] = sinthetaj * cos(phik);
            gridPoints[3 * index + 1] = sinthetaj * sin(phik);
            gridPoints[3 * index + 2] = costhetaj;
            gridWeights[index] = weightfactor * weightsGL[p - j]; // area element = sin thetaj
        }
    }
    // south pole
    gridPoints[3 * (Ngrid - 1)] = 0;
    gridPoints[3 * (Ngrid - 1) + 1] = 0;
    gridPoints[3 * (Ngrid - 1) + 2] = -1;
    gridWeights[Ngrid - 1] = 0;

    // rotation with quaternion for each point coordinate
    for (int i = 0; i < Ngrid; i++) {
        Evec3 pointVec(gridPoints[3 * i], gridPoints[3 * i + 1], gridPoints[3 * i + 2]);
        Eigen::Map<Evec3>(gridPoints.data() + 3 * i) = (orientation * pointVec);
    }

    // fill grid norms
    if (gridNormPtr != nullptr) {
        *gridNormPtr = gridPoints; // on unit sphere grid norms equal grid coordinates
    }

    // scale and shift grid points
    for (int i = 0; i < Ngrid; i++) {
        Evec3 pointVec(gridPoints[3 * i], gridPoints[3 * i + 1], gridPoints[3 * i + 2]);
        Eigen::Map<Evec3>(gridPoints.data() + 3 * i) = (radius * pointVec) + coordBase;
    }

    return;
}

void Shexp::getGridWithPoleCellConnect(std::vector<int32_t> &gridCellConnect, std::vector<int32_t> &offset,
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

// excluding the poles
// if valPtr=nullptr, directly fill the gridValues
// spectralCoeff pointed by coeffPtr is not modified
void Shexp::calcGridValue(double *coeffPtr, double *valPtr) {

    if (valPtr == nullptr) {
        valPtr = gridValue.data();
    }

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
        // rearrange the grid data from  (x0,x1,...y0,y1,...z0,z1...) to (x0,y0,z0,x1,y1,z1...)
        std::vector<double> gridValueBuffer(Ngrid * dimension);
        for (int i = 0; i < Ngrid; i++) {
            gridValueBuffer[3 * i] = valPtr[i];                 // xi
            gridValueBuffer[3 * i + 1] = valPtr[Ngrid + i];     // yi
            gridValueBuffer[3 * i + 2] = valPtr[2 * Ngrid + i]; // zi
        }
        rotGridValue(gridValueBuffer.data(), Ngrid);
        for (int i = 0; i < 3 * Ngrid; i++) {
            valPtr[i] = gridValueBuffer[i];
        }
    }
}

// spectralCoeff pointed by coeffPtr is not modified
void Shexp::calcPoleValue(double *coeffPtr, double *valPtr) const {

    typedef double Real;
    const int Ncoeff = (order + 1) * (order + 2);
    const int Ngrid = 2; // north and south pole
    const int dimension = getDimension();

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> val(Ngrid * dimension, sctl::Ptr2Itr<Real>(valPtr ? valPtr : nullptr, Ngrid * dimension), false);

    sctl::Vector<Real> cos_theta_phi(4);
    cos_theta_phi[0] = 1.0;  // north pole, cos theta
    cos_theta_phi[1] = 0;    // north pole, phi
    cos_theta_phi[2] = -1.0; // south pole, cos theta
    cos_theta_phi[3] = 0;    // south pole, phi

    if (kind == KIND::LAP) {
        sctl::SphericalHarmonics<Real>::SHCEval(coeff, sctl::SHCArrange::ROW_MAJOR, order, cos_theta_phi, val);
    } else {
        sctl::SphericalHarmonics<Real>::VecSHCEval(coeff, sctl::SHCArrange::ROW_MAJOR, order, cos_theta_phi, val);
        rotGridValue(valPtr, 2);
    }
}

// gridValue pointed by valPtr is not modified
void Shexp::calcSpectralCoeff(double *coeffPtr, const double *valPtr) const {

    typedef double Real;
    const int Ncoeff = (order + 1) * (order + 2);
    const int Ngrid = (order + 1) * (2 * order + 2);
    const int dimension = getDimension();

    sctl::Vector<Real> coeff(Ncoeff * dimension, sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension),
                             false);
    sctl::Vector<Real> val(Ngrid * dimension);

    // at this point, val and valPtr always have same values

    if (kind == KIND::LAP) {
        if (valPtr == nullptr) {
            val = gridValue;
        } else {
            std::copy(valPtr, valPtr + Ngrid, val.begin());
        }
        sctl::SphericalHarmonics<Real>::Grid2SHC(val, order + 1, 2 * order + 2, order, coeff,
                                                 sctl::SHCArrange::ROW_MAJOR);
    } else {
        std::vector<double> gridValueBuffer; // temporary space for rotation
        if (valPtr == nullptr) {
            gridValueBuffer = gridValue;
        } else {
            gridValueBuffer.resize(Ngrid * dimension);
            std::copy(valPtr, valPtr + Ngrid * dimension, gridValueBuffer.begin());
        }
        invrotGridValue(gridValueBuffer.data(), Ngrid);
        // rearrange the grid data from (x0,y0,z0,x1,y1,z1...) to  (x0,x1,...y0,y1,...z0,z1...)
        for (int i = 0; i < Ngrid; i++) {
            val[i] = gridValueBuffer[3 * i];
            val[Ngrid + i] = gridValueBuffer[3 * i + 1];
            val[2 * Ngrid + i] = gridValueBuffer[3 * i + 2];
        }
        sctl::SphericalHarmonics<Real>::Grid2VecSHC(val, order + 1, 2 * order + 2, order, coeff,
                                                    sctl::SHCArrange::ROW_MAJOR);
    }
}

void Shexp::rotGridValue(double *valPtr, const int npts) const {
    if (orientation * Evec3(0, 0, 1) == Evec3(0, 0, 1)) {
        return;
    }

    // from default Z axis to oriented Z axis
    assert(kind == KIND::STK);
    // rotation with quaternion for each point vectorial value
    Equatn q = orientation;

    for (int i = 0; i < npts; i++) {
        Evec3 vec(valPtr[3 * i], valPtr[3 * i + 1], valPtr[3 * i + 2]);
        Emap3(valPtr + 3 * i) = q * vec;
    }

    return;
}
void Shexp::invrotGridValue(double *valPtr, const int npts) const {
    if (orientation * Evec3(0, 0, 1) == Evec3(0, 0, 1)) {
        return;
    }
    // from oriented Z axis to default Z axis
    assert(kind == KIND::STK);
    // rotation with quaternion for each point vectorial value
    Equatn q = orientation.inverse();

    for (int i = 0; i < npts; i++) {
        Evec3 vec(valPtr[3 * i], valPtr[3 * i + 1], valPtr[3 * i + 2]);
        Emap3(valPtr + 3 * i) = q * vec;
    }
    return;
}

void Shexp::randomFill(const int seed) {
    std::vector<double> spectralCoeff(getSpectralDOF(), 0);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-2.0, 2.0);

    const int Ncoeff = (order + 1) * (order + 2); // Ncoeff doubles, represent Ncoeff/2 complex

    for (int d = 0; d < getDimension(); d++) {
        for (int i = 0; i < Ncoeff; i++) {
            spectralCoeff[i + d * Ncoeff] = dis(gen) * pow(0.8, i);
        }
    }

    calcGridValue(spectralCoeff.data());
    // 1 more transform to remove complex components
    calcSpectralCoeff(spectralCoeff.data());
    calcGridValue(spectralCoeff.data());

    return;
}

void Shexp::calcSDLNF(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgValuePtr, const bool &interior,
                      const bool &SL) const {
    using Real = double;
    const int Ncoeff = (order + 1) * (order + 2);
    const int dimension = getDimension();

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> xyz(trgNum * 3, sctl::Ptr2Itr<Real>(trgXYZPtr ? trgXYZPtr : nullptr, trgNum * 3), false);

    sctl::Vector<Real> val(trgNum * dimension,
                           sctl::Ptr2Itr<Real>(trgValuePtr ? trgValuePtr : nullptr, trgNum * dimension), false);

    const double invRadius = 1.0 / radius;
    xyz *= invRadius; // scale with radius

    // rotate XYZ Coord to sh's frame;
    invrotGridValue(trgXYZPtr, trgNum);

    // compute & output
    if (dimension == 1) {
        // no rotation of grid value needed
        // compute
        printf("not implemented\n");
        exit(1);
    } else {
        // compute
        if (SL) {
            sctl::SphericalHarmonics<Real>::StokesEvalSL(coeff, sctl::SHCArrange::ROW_MAJOR, order, xyz, interior, val);
            val *= (radius); // scale back
        } else {
            sctl::SphericalHarmonics<Real>::StokesEvalDL(coeff, sctl::SHCArrange::ROW_MAJOR, order, xyz, interior, val);
            // radius scaling already satisfied
            val *= -1; // prefactor set to -3/4pi
        }
        // rotate back to lab frame for stk
        rotGridValue(trgValuePtr, trgNum);
    }

    return;
}

void Shexp::calcKSelf(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgValuePtr,
                      const bool &interior) const {
    using Real = double;
    const int Ncoeff = (order + 1) * (order + 2);
    const int dimension = getDimension();

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> xyz(trgNum * 3, sctl::Ptr2Itr<Real>(trgXYZPtr ? trgXYZPtr : nullptr, trgNum * 3), false);

    sctl::Vector<Real> val(trgNum * dimension,
                           sctl::Ptr2Itr<Real>(trgValuePtr ? trgValuePtr : nullptr, trgNum * dimension), false);

    const double invRadius = 1.0 / radius;
    xyz *= invRadius; // scale with radius

    // rotate XYZ Coord to sh's frame;
    invrotGridValue(trgXYZPtr, trgNum);

    // compute & output
    if (dimension == 1) {
        // no rotation of grid value needed
        // compute
        printf("error, K is for Stokes only\n");
        exit(1);
    } else {
        // compute
        sctl::SphericalHarmonics<Real>::StokesEvalKSelf(coeff, sctl::SHCArrange::ROW_MAJOR, order, xyz, interior, val);
        // radius scaling already satisfied
        // rotate back to lab frame for stk
        rotGridValue(trgValuePtr, trgNum);
    }

    return;
}

void Shexp::calcKNF(double *coeffPtr, const int &trgNum, double *trgXYZPtr, double *trgNormPtr, double *trgValuePtr,
                    const bool &interior) const {
    using Real = double;
    const int Ncoeff = (order + 1) * (order + 2);
    const int dimension = 3;

    const sctl::Vector<Real> coeff(Ncoeff * dimension,
                                   sctl::Ptr2Itr<Real>(coeffPtr ? coeffPtr : nullptr, Ncoeff * dimension), false);
    sctl::Vector<Real> xyz(trgNum * 3, sctl::Ptr2Itr<Real>(trgXYZPtr ? trgXYZPtr : nullptr, trgNum * 3), false);
    sctl::Vector<Real> norm(trgNum * 3, sctl::Ptr2Itr<Real>(trgNormPtr ? trgNormPtr : nullptr, trgNum * 3), false);

    sctl::Vector<Real> val(trgNum * dimension,
                           sctl::Ptr2Itr<Real>(trgValuePtr ? trgValuePtr : nullptr, trgNum * dimension), false);

    const double invRadius = 1.0 / radius;
    xyz *= invRadius; // scale with radius

    // rotate XYZ Coord to sh's frame;
    invrotGridValue(trgXYZPtr, trgNum);
    invrotGridValue(trgNormPtr, trgNum);

    // printf("%lf,%lf,%lf\n", trgXYZPtr[0], trgXYZPtr[1], trgXYZPtr[2]);
    // printf("%lf,%lf,%lf\n", trgNormPtr[0], trgNormPtr[1], trgNormPtr[2]);

    // compute & output
    if (dimension == 1) {
        // no rotation of grid value needed
        // comput
        printf("not implemented\n");
        exit(1);
    } else {
        // compute
        sctl::SphericalHarmonics<Real>::StokesEvalKL(coeff, sctl::SHCArrange::ROW_MAJOR, order, xyz, norm, interior,
                                                     val);
        // printf("%lf,%lf,%lf\n", trgValuePtr[0], trgValuePtr[1], trgValuePtr[2]);

        // radius scaling already satisfied
        // rotate back to lab frame for stk
        // invrotGridValue(trgValuePtr, trgNum);
        rotGridValue(trgValuePtr, trgNum);
    }

    return;
}