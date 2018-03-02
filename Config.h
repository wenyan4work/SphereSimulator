#ifndef CONFIG_H
#define CONFIG_H

#include <string>

class Config {
  public:
    // parallel setting
    int ompThreads;

    // domain setting
    double simBoxHigh[3];                // simulation box size
    double simBoxLow[3];                 // simulation box size
    int xPeriodic, yPeriodic, zPeriodic; // flag of true/false of periodic in that direction
    bool monolayer;
    bool hydro;
    double scaleBrown;
    double StkReg;
    bool dumpflow;
    bool shell;
    int pFMM;            // mult_order p for PVFMM
    double dumpFlowMesh; // flow dump mesh size

    // physical setting
    double sphereRadiusHydro;          // um
    double sphereRadiusSigmaHydro;     // sigma for log normal distribution
    double sphereRadiusCollisionRatio; // ratio*radiusHydro=radiusCollision

    // physical constant
    double viscosity;       // pN/(um^2 s)
    double kBT;             // pN.um

    // number
    int sphereNumber;
    unsigned int rngSeed;

    // time stepping
    double dt;
    double timeTotal;
    int snapFreq;

    explicit Config(std::string);
    ~Config() = default;
};

#endif