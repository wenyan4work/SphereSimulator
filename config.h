#include <string>

class Config {
public:
    // parallel setting
    int xProcs, yProcs, zProcs; // number of processes in each direction
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
    int pFMM; // mult_order p for PVFMM
    double dumpFlowMesh; // flow dump mesh size

    // physical setting
    double sphereRadiusA;      // um
    double sphereRadiusSigmaA; // sigma for log normal distribution
    double sphereRadiusB;      // um
    double sphereRadiusSigmaB; // sigma for log normal distribution

    // physical constant
    double viscosity;       // pN/(um^2 s)
    double kBT;             // pN.um
    double swimForce;       // pN
    double springLength;    // um
    double springMaxLength; // um
    double forceMaxRatio;   // the ratio of max force to swim force

    // number
    int sphereNumber;
    unsigned int rngSeed;

    // time stepping
    double dt;
    double timeTotal;
    double timeSnap;

    explicit Config(std::string);
    ~Config() = default;
};
