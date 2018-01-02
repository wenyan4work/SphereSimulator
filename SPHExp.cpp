#include <cstdio>

#include "SPHExp.hpp"

// constructor
SPHExp::SPHExp(const KIND kind_, const std::string &name_, const int order_, const Equatn orientation_)
    : kind(kind_), order(order_), dimension((kind_ == LAPDL || kind_ == LAPSL) ? 1 : 3), name(name_),
      orientation(orientation_) { // the dimension , =1 for LAPSL, LAPDL, =3 for STKSL and STKDL

    spectralCoeff.resize(dimension * (order + 1) * (order + 1)); // coefficients. d*(p+1)^2 elements
}

// copy constructor
SPHExp::SPHExp(const SPHExp &other)
    : kind(other.kind), order(other.order), dimension((other.kind == LAPDL || other.kind == LAPSL) ? 1 : 3),
      name(other.name), orientation(other.orientation) {
    spectralCoeff = other.spectralCoeff;
}

SPHExp::SPHExp(SPHExp &&other)
    : kind(other.kind), order(other.order), dimension((other.kind == LAPDL || other.kind == LAPSL) ? 1 : 3),
      name(other.name), orientation(other.orientation) {
    spectralCoeff = std::move(other.spectralCoeff);
}

void SPHExp::dumpSpectralValues(const std::string &filename) {
    std::string kindName;
    switch (kind) {
    case LAPSL:
        kindName = "LapSL";
        break;
    case LAPDL:
        kindName = "LapDL";
        break;
    case STKSL:
        kindName = "StkSL";
        break;
    case STKDL:
        kindName = "StkDL";
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
