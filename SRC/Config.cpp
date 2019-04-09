#include "Config.h"

#include <fstream>
#include <iostream>
#include <omp.h>

void readline(FILE *file) {
    while (fgetc(file) != '\n') {
    }
    return;
}

Config::Config(std::string fileName) {
    std::ifstream myfile(fileName);
    if (!myfile.is_open()) {
        printf("Error: cannot open runConfig.txt file");
        // MPI_Finalize();
        exit(1);
    }
    std::string line;

    getline(myfile, line);
    ompThreads = 1;
    myfile >> ompThreads >> line;
    if (ompThreads > 0) {
        omp_set_num_threads(ompThreads);
        omp_set_nested(0);
    }
    std::cout << "using max OpenMP threads: " << omp_get_max_threads() << std::endl;

    getline(myfile, line);
    simBoxLow[0] = 0;
    simBoxLow[1] = 0;
    simBoxLow[2] = 0;
    simBoxHigh[0] = 0;
    simBoxHigh[1] = 0;
    simBoxHigh[2] = 0;
    myfile >> simBoxLow[0] >> simBoxLow[1] >> simBoxLow[2] >> simBoxHigh[0] >> simBoxHigh[1] >> simBoxHigh[2] >> line;

    getline(myfile, line);
    xPeriodic = 0;
    yPeriodic = 0;
    zPeriodic = 0;
    myfile >> xPeriodic >> yPeriodic >> zPeriodic >> line;

    getline(myfile, line);
    sphereRadiusHydro = 0;
    sphereRadiusSigmaHydro = 0;
    myfile >> sphereRadiusHydro >> sphereRadiusSigmaHydro >> line;

    getline(myfile, line);
    sphereRadiusCollisionRatio = 0;
    myfile >> sphereRadiusCollisionRatio >> line;

    getline(myfile, line);
    viscosity = 0;
    myfile >> viscosity >> line;

    getline(myfile, line);
    kBT = 0;
    myfile >> kBT >> line;

    getline(myfile, line);
    extForce[0] = extForce[1] = extForce[2] = 0;
    extTorque[0] = extTorque[1] = extTorque[2] = 0;
    myfile >> extForce[0] >> extForce[1] >> extForce[2] >> extTorque[0] >> extTorque[1] >> extTorque[2] >> line;

    getline(myfile, line);
    sphereNumber = 0;
    myfile >> sphereNumber >> line;

    getline(myfile, line);
    dt = 0;
    myfile >> dt >> line;

    getline(myfile, line);
    timeTotal = 0;
    myfile >> timeTotal >> line;

    getline(myfile, line);
    snapFreq = 0;
    myfile >> snapFreq >> line;

    getline(myfile, line);
    rngSeed = 0;
    myfile >> rngSeed >> line;

    getline(myfile, line);
    monolayer = true;
    int temp = 0;
    myfile >> temp >> line;
    monolayer = temp > 0 ? true : false;

    getline(myfile, line);
    hydro = true;
    temp = 0;
    myfile >> temp >> line;
    hydro = temp > 0 ? true : false;

    getline(myfile, line);
    pSH = 3;
    std::cout << line;
    myfile >> pSH >> line;

    getline(myfile, line);
    pFMM = 6;
    std::cout << line;
    myfile >> pFMM >> line;
    if (pFMM % 2 != 0) {
        printf("P for FMM must be an even number: p=%d\n", pFMM);
        exit(1);
    } else if (pFMM < 6 || pFMM > 16) {
        printf("P for FMM must be between [6,16]: p=%d\n", pFMM);
        exit(1);
    } else {
        printf("using pFMM %d\n", pFMM);
    }

    getline(myfile, line);
    scaleBrown = 1.0;
    double tempd = 0;
    myfile >> tempd >> line;
    if (tempd < 0) {
        printf("scale Brown must be nonnegative. using scaleBrown = 1.");
    } else {
        scaleBrown = tempd;
        printf("using scaleBrown = %lf.", scaleBrown);
    }

    getline(myfile, line);
    StkReg = 0;
    myfile >> StkReg >> line;
    if (StkReg > 0) {
        std::cout << "Using Stokes Regularization Parameter " << StkReg << std::endl;
    }

    getline(myfile, line);
    dumpflow = false;
    myfile >> temp >> line;
    if (temp > 0) {
        std::cout << "generate flow maps" << std::endl;
        dumpflow = true;
    }

    getline(myfile, line);
    dumpFlowMesh = 0;
    myfile >> dumpFlowMesh >> line;
    if (dumpFlowMesh < 0) {
        std::cout << "Dump flow mesh size must be positive" << std::endl;
        exit(1);
    }

    getline(myfile, line);
    shell = false;
    myfile >> temp >> line;
    if (temp > 0) {
        std::cout << "use spherical shell as boundary" << std::endl;
        shell = true;
    }

    myfile.close();

    // input correctness check
    {
        printf("Run Setting: \n");
        printf("Simulation box Low: %lf,%lf,%lf\n", simBoxLow[0], simBoxLow[1], simBoxLow[2]);
        printf("Simulation box High: %lf,%lf,%lf\n", simBoxHigh[0], simBoxHigh[1], simBoxHigh[2]);
        printf("Periodicity: %d,%d,%d\n", xPeriodic > 0 ? true : false, yPeriodic > 0 ? true : false,
               zPeriodic > 0 ? true : false);
        printf("Monolayer: %d\n", monolayer);
        printf("With hydrodynamics: %d\n", hydro);
    }
    {
        printf("Sphere Setting: \n");
        printf("Sphere radius average: %lf, LogNormal sigma: %lf \n", sphereRadiusHydro, sphereRadiusSigmaHydro);
        printf("Sphere radius collision ratio: %lf\n", sphereRadiusCollisionRatio);
    }

    {
        printf("Physical setting: \n");
        printf("viscosity: %lf\n", viscosity);
        printf("kBT: %lf\n", kBT);
    }
    {
        printf("Sphere number: %d\n", sphereNumber);
        printf("Time step size: %lf\n", dt);
        printf("Total Time: %lf\n", timeTotal);
        printf("Snap Freq: %d\n", snapFreq);
    }

    return;
}
