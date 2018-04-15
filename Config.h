#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include "Util/EigenDef.hpp"

/**
 * \class Config
 *
 * \brief Set up the simulation configuration
 */
class Config {
  public:
    /// How many OpenMP threads will be used per process
    int ompThreads;

    /// Upper corner of the simulation box
    double simBoxHigh[3];
    /// Lower corner of the simulation box
    double simBoxLow[3];
    /// Whether x direction is periodic
    int xPeriodic;
    /// Whether y direction is periodic
    int yPeriodic;
    /// Whether z direction is periodic
    int zPeriodic;
    ///
    bool monolayer;
    /// Whether fluid layers are enabled
    bool hydro;
    /// 
    double scaleBrown;
    ///
    double StkReg;
    ///
    bool dumpflow;
    /// Whether a spherical shell is present
    bool shell;
    /// Shell center
    Evec3 shellCenter;
    /// Shell radius
    double shellRadius;
    /// Mutlipole order p for PVFMM
    int pFMM;
    /// flow dump mesh size (?)
    double dumpFlowMesh;

    /// Radius of sphere (um)
    double sphereRadiusHydro;
    /// Variance for log-normal distribution to generate random radius
    double sphereRadiusSigmaHydro;
    /// Ratio of collision radius and sphere radius
    double sphereRadiusCollisionRatio;

    /// Acceleration due to gravity (um/s)
    double gravity;
    /// Coeffecient of Viscosity (pN/um^2-s)
    double viscosity;
    /// Product of Boltzmann constant and temperature (pN-um)
    double kBT;

    /// Total number of spheres
    int sphereNumber;
    /// Random number seed
    unsigned int rngSeed;

    /// Timestep
    double dt;
    /// Total simulation time
    double timeTotal;
    /// How often to write simulation results to disk
    int snapFreq;

    /**
     * \brief Read simulation configuration from file
     */
    explicit Config(std::string);
    ~Config() = default;
};

#endif
