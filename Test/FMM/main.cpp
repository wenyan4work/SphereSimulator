/*
 * main.cpp
 *
 *  Created on: Oct 14, 2016
 *      Author: wyan
 */

#include "cmdparser.hpp"
#include "fmm.h"
#include "mpi.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

extern "C" {
void stfmm3dparttarg_(int *ier, int *iprec, int *nsrc, double src[], int *ifsingle, double sigma_sl[], int *ifdouble,
                     double sigma_dl[], double sigma_dv[], int *ifpot, double pot[], double pre[], int *ifgrad,
                     double grad[], int *ntrg, double trg[], int *ifpottrg, double pottrg[], double pretrg[],
                     int *ifgradtrg, double gradtrg[]);
}

void calcSTFMM(std::vector<double> &src_coord, std::vector<double> &srcSL, std::vector<double> &srcDL,
               std::vector<double> &srcDV, std::vector<double> &srcPre, std::vector<double> &srcPot,
               std::vector<double> &srcPotGrad, std::vector<double> &trg_coord, std::vector<double> &trgPre,
               std::vector<double> &trgPot, std::vector<double> &trgPotGrad) {
    // call the library STFMM3D
    int ier = 0;
    int iprec = 1;
    int nsrc = src_coord.size() / 3;
    int ifsingle = 1;
    int ifdouble = 1;
    int ifpot = 1;
    int ifgrad = 1;
    int ntrg = trg_coord.size() / 3;
    int ifpottrg = 1;
    int ifgradtrg = 1;

    stfmm3dparttarg_(&ier, &iprec, &nsrc, src_coord.data(), &ifsingle, srcSL.data(), &ifdouble, srcDL.data(),
                    srcDV.data(), &ifpot, srcPot.data(), srcPre.data(), &ifgrad, srcPotGrad.data(), &ntrg,
                    trg_coord.data(), &ifpottrg, trgPot.data(), trgPre.data(), &ifgradtrg, trgPotGrad.data());

    for (auto &v : trgPre) {
        std::cout << "trgPre" << v << std::endl;
    }
    for (auto &v : trgPot) {
        std::cout << "trgPot" << v << std::endl;
    }
    for (auto &v : trgPotGrad) {
        std::cout << "trgPotGrad" << v << std::endl;
    }

    return;
}

void calcPVFMM(std::vector<double> &src_coord, std::vector<double> &srcSL, std::vector<double> &srcDL,
               std::vector<double> &srcDV, std::vector<double> &srcPre, std::vector<double> &srcPot,
               std::vector<double> &srcPotGrad, std::vector<double> &trg_coord, std::vector<double> &trgPre,
               std::vector<double> &trgPot, std::vector<double> &trgPotGrad) {
    // call the library STFMM3D
    int ier = 0;
    int iprec = 1;
    int nsrc = src_coord.size() / 3;
    int ifsingle = 1;
    int ifdouble = 1;
    int ifpot = 1;
    int ifgrad = 1;
    int ntrg = trg_coord.size() / 3;
    int ifpottrg = 1;
    int ifgradtrg = 1;

    stfmm3dparttarg_pvfmm(&ier, &iprec, &nsrc, src_coord.data(), &ifsingle, srcSL.data(), &ifdouble, srcDL.data(),
                          srcDV.data(), &ifpot, srcPot.data(), srcPre.data(), &ifgrad, srcPotGrad.data(), &ntrg,
                          trg_coord.data(), &ifpottrg, trgPot.data(), trgPre.data(), &ifgradtrg, trgPotGrad.data());

    for (auto &v : trgPre) {
        std::cout << "pvfmm trgPre" << v << std::endl;
    }
    for (auto &v : trgPot) {
        std::cout << "pvfmm trgPot" << v << std::endl;
    }
    for (auto &v : trgPotGrad) {
        std::cout << "pvfmm trgPotGrad" << v << std::endl;
    }

    return;
}

void compare(const std::vector<double> &trg_value, const std::vector<double> &trg_value_true) {
    std::cout << "result, fmm, true:" << std::endl;
    for (int i = 0; i < trg_value_true.size(); i++) {
        std::cout << "pvfmm " << trg_value[i] << " stfmm3d " << trg_value_true[i] << std::endl;
    }

    // calc error and max error
    double errorL2 = 0, errorAbs = 0, L2 = 0, errorMaxL2 = 0, maxU = 0;
    double errorMaxRel = 0;
    for (int i = 0; i < trg_value_true.size(); i++) {
        double temp = pow(trg_value_true[i] - trg_value[i], 2);
        //		if (temp >= pow(1e-5, 2)) {
        //			std::cout << "i" << i << "error L2" << temp <<
        // std::endl;
        //		}
        errorL2 += temp;
        errorAbs += sqrt(temp);
        L2 += pow(trg_value_true[i], 2);
        errorMaxL2 = std::max(sqrt(temp), errorMaxL2);
        maxU = std::max(maxU, fabs(trg_value_true[i]));
        errorMaxRel = std::max(sqrt(temp) / trg_value_true[i], errorMaxRel);
    }

    std::cout << std::setprecision(16) << "Max Abs Error L2: " << (errorMaxL2) << std::endl;
    std::cout << std::setprecision(16) << "Ave Abs Error L2: " << errorAbs / trg_value_true.size() << std::endl;
    std::cout << std::setprecision(16) << "Max Rel Error L2: " << (errorMaxRel) << std::endl;
    std::cout << std::setprecision(16) << "RMS Error L2: " << sqrt(errorL2 / trg_value_true.size()) << std::endl;
    std::cout << std::setprecision(16) << "Relative Error L2: " << sqrt(errorL2 / L2) << std::endl;
}

void configure_parser(cli::Parser &parser) {
    parser.set_optional<int>("T", "ntarget", 2, "target number in each dimension. default 2");
    parser.set_optional<double>("B", "box", 1.0, "box edge length");
    parser.set_optional<double>("M", "move", 0.0, "box origin shift move");
    parser.set_optional<int>("S", "source", 1, "1 for point force, 2 for force dipole, other for same as target.");
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    cli::Parser parser(argc, argv);
    configure_parser(parser);
    parser.run_and_exit_if_error();

    std::cout << "Running setting: " << std::endl;
    std::cout << "nTarget: " << parser.get<int>("T") << std::endl;
    std::cout << "Box: " << parser.get<double>("B") << std::endl;
    std::cout << "Shift: " << parser.get<double>("M") << std::endl;
    std::cout << "Source: " << parser.get<int>("S") << std::endl;

    // FMM_Wrapper::PAXIS pset = FMM_Wrapper::PAXIS::NONE;

    const double box = parser.get<double>("B");
    const double shift = parser.get<double>("M");
    const int ntrg = parser.get<int>("T");

    std::vector<double> src_coord(0);
    std::vector<double> trg_coord(0);

    // initialize source and target coord and value

    // initialize target points
    std::random_device rd;
    std::mt19937 gen(rd());

    std::lognormal_distribution<> d(log(0.2), 0.5);

    for (int i = 0; i < ntrg; i++) {
        double r01 = fmod(d(gen), 1.0);
        trg_coord.push_back(r01 * box + shift); // x
        r01 = fmod(d(gen), 1.0);
        trg_coord.push_back(r01 * box + shift); // y
        r01 = fmod(d(gen), 1.0);
        trg_coord.push_back(r01 * box + shift); // z
    }

    FILE *pfile = fopen("randomPoints", "w");
    for (int i = 0; i < trg_coord.size() / 3; i++) {
        fprintf(pfile, "%f\t%f\t%f\n", trg_coord[3 * i], trg_coord[3 * i + 1], trg_coord[3 * i + 2]);
    }
    fclose(pfile);

    // initialize source points
    int nsrc = parser.get<int>("S");
    for (int i = 0; i < nsrc; i++) {
        double r01 = fmod(d(gen), 1.0);
        src_coord.push_back(r01 * box + shift); // x
        r01 = fmod(d(gen), 1.0);
        src_coord.push_back(r01 * box + shift); // y
        r01 = fmod(d(gen), 1.0);
        src_coord.push_back(r01 * box + shift); // z
    }

    // check size
    const int n_src = src_coord.size() / 3;
    const int n_trg = trg_coord.size() / 3;
    if (n_src * 3 != src_coord.size()) {
        std::cout << "src size error" << std::endl;
        exit(1);
    }
    if (n_trg * 3 != trg_coord.size()) {
        std::cout << "trg size error" << std::endl;
        exit(1);
    }

    for (auto &s : src_coord) {
        std::cout << "src_coord" << s << std::endl;
    }
    for (auto &t : trg_coord) {
        std::cout << "trg_coord" << t << std::endl;
    }
    // setup source
    std::vector<double> srcSL(n_src * 3); // single layer
    std::vector<double> srcDL(n_src * 3); // double layer
    std::vector<double> srcDV(n_src * 3); // double layer orientation. no need to normalize
    for (auto &v : srcSL) {
        v = d(gen);
        v = 0;
        std::cout << "srcSL" << v << std::endl;
    }
    for (auto &v : srcDL) {
        v = d(gen);
        std::cout << "srcDL" << v << std::endl;
    }
    for (auto &v : srcDV) {
        v = d(gen);
        std::cout << "srcDV" << v << std::endl;
    }

    std::vector<double> srcPre(n_src);
    std::vector<double> srcPot(n_src * 3);
    std::vector<double> srcPotGrad(n_src * 9);
    std::fill(srcPre.begin(), srcPre.end(), 0);
    std::fill(srcPot.begin(), srcPot.end(), 0);
    std::fill(srcPotGrad.begin(), srcPotGrad.end(), 0);
    std::vector<double> trgPre(n_trg);
    std::vector<double> trgPot(n_trg * 3);
    std::vector<double> trgPotGrad(n_trg * 9);
    std::fill(trgPre.begin(), trgPre.end(), 0);
    std::fill(trgPot.begin(), trgPot.end(), 0);
    std::fill(trgPotGrad.begin(), trgPotGrad.end(), 0);

    // calc true value
    std::vector<double> srcPreTrue = srcPre;
    std::vector<double> srcPotTrue = srcPot;
    std::vector<double> srcPotGradTrue = srcPotGrad;

    std::vector<double> trgPreTrue = trgPre;
    std::vector<double> trgPotTrue = trgPot;
    std::vector<double> trgPotGradTrue = trgPotGrad;

    calcSTFMM(src_coord, srcSL, srcDL, srcDV, srcPreTrue, srcPotTrue, srcPotGradTrue, trg_coord, trgPreTrue, trgPotTrue,
              trgPotGradTrue);

    calcPVFMM(src_coord, srcSL, srcDL, srcDV, srcPre, srcPot, srcPotGrad, trg_coord, trgPre, trgPot, trgPotGrad);

    printf("srcPre\n");
    compare(srcPre, srcPreTrue);

    printf("srcPot\n");
    compare(srcPot, srcPotTrue);

    printf("srcPotGrad\n");
    compare(srcPotGrad, srcPotGradTrue);

    printf("trgPre\n");
    compare(trgPre, trgPreTrue);

    printf("trgPot\n");
    compare(trgPot, trgPotTrue);

    printf("trgPotGrad\n");
    compare(trgPotGrad, trgPotGradTrue);

    MPI_Finalize();

    return 0;
}
