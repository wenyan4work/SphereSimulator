//#pragma GCC push_options
//// for debug
//#pragma GCC optimize ("O0")
#ifndef SPHERESYSTEM_H_
#define SPHERESYSTEM_H_

#include "FDPS_Header/particle_simulator.hpp"
#include "TRngPool.hpp"
#include "config.h"
#include "LayerSphereSystem.hpp"

#include <algorithm>
#include <deque>
#include <memory>
#include <set>

#define COLBUF 0.3 // possible collision within (1+0.3)*(radiusI+radiusJ)

// check whether an integer is even or odd
inline bool isOdd(const int x) { return x % 2 != 0; }

inline double normVec3(const PS::F64vec3 &vec3) {
    double norm2 = vec3 * vec3;
    return sqrt(norm2);
}

inline void normalizeVec3(PS::F64vec3 &vec3) {
    double norm2 = vec3 * vec3;
    double rnorm = 1 / sqrt(norm2);
    vec3.x *= rnorm;
    vec3.y *= rnorm;
    vec3.z *= rnorm;
}

inline void zeroVec3(PS::F64vec3 &vec3) {
    vec3.x = 0;
    vec3.y = 0;
    vec3.z = 0;
}

// output file header
class FileHeader {
  public:
    int nSphere;
    PS::F64 time;

    void writeAscii(FILE *fp) const {
        fprintf(fp, "%d\n", nSphere);
        fprintf(fp, "SphereNumber: %d, time: %lf\n", nSphere, time);
    }
};

/*FDPS force class*/

// type for near pairwise force. not useful in collision detection
class forceNearEP {
  public:
    PS::F64vec3 forceNear;
    PS::F64vec3 torqueNear;
    PS::F64vec3 sepmin;

    void clear() {
        forceNear = 0;
        torqueNear = 0;
        sepmin = 0;
    }
};

/*FullParticle class for FDPS*/
class SphereFP {
  public:
    PS::S32 gid, localIndex, globalIndex;
    // PS::S32 rank;
    PS::F64vec3 pos; // Tubule center or protein center
    PS::F64 radius;
    PS::F64 RSearch;
    // the search radius for gather and scatter operations, the actually used radius are in the EP classes
    // FDPS support different search radius for different TYPE of EPI/EPJ,
    // but does NOT support different search radius for different ENTITY of particle.
    // radius Search should be set to the LONGEST one for all entities of the same type.
    // line 881 of FDPS_tree_for_force_impl.hpp:
    // const F64 r_cut_sq  = epj_sorted_[0].getRSearch() * epj_sorted_[0].getRSearch();

    PS::F64vec3 randTrans; // unit rng N01
    PS::F64vec3 randRot;   // unit rng N01
    PS::F64vec3 posT0;

    PS::F64vec3 forceExt; // external
    PS::F64vec3 torqueExt;

    PS::F64vec3 forceNear; // Near
    PS::F64vec3 torqueNear;

    PS::F64vec3 forceCollision; // collision resolving
    PS::F64vec3 torqueCollision;

    PS::F64vec3 vel;   // velocity
    PS::F64vec3 omega; // angular velocity

    PS::S32 overlap;

    PS::F64vec getPos() const { return this->pos; }

    void copyFromForce(const forceNearEP &force) {
        this->forceNear = force.forceNear;
        this->torqueNear = force.torqueNear;
    }

    void setPos(const PS::F64vec3 &pos_) {
        this->posT0 += (pos_ - this->pos); // posT0 is the 'relative' initial position
        this->pos = pos_;
    }

    void writeAscii(FILE *fp) const {
        fprintf(fp, "S\t%8d\t%.17g\t%.17g\t%.17g\t%.17g\n", this->gid, this->radius, this->pos[0], this->pos[1],
                this->pos[2]);
    }
};

class SphereNearForceEP {
  public:
    PS::S32 gid, localIndex, globalIndex;

    // PS::S32 rank; // on which rank

    PS::F64vec3 pos;
    PS::F64 radius;
    PS::F64 RSearch;

    PS::F64vec getPos() const { return this->pos; }

    void setPos(const PS::F64vec3 &pos) { this->pos = pos; }

    PS::F64 getRSearch() const { return this->RSearch; }

    void copyFromFP(const SphereFP &fp) {
        gid = fp.gid;
        localIndex = fp.localIndex;
        globalIndex = fp.globalIndex;
        pos = fp.pos;
        radius = fp.radius;
        RSearch = fp.RSearch;
    }
};

struct ColBlock {
  public:
    int gidI, gidJ;
    // int localIndexI, localIndexJ;
    int globalIndexI, globalIndexJ;
    PS::F64vec3 gvecI, gvecJ;
    PS::F64vec3 posI, posJ;
    double phi0;   // constraint value
    double lambda; // force magnitude

    ColBlock() : gidI(0), gidJ(0), globalIndexI(0), globalIndexJ(0), phi0(0), lambda(0) {
        // default constructor
        zeroVec3(gvecI);
        zeroVec3(gvecJ);
        zeroVec3(posI);
        zeroVec3(posJ);
    }

    ColBlock(int gidI_, int gidJ_, int globalIndexI_, int globalIndexJ_, const PS::F64vec3 &gvecI_,
             const PS::F64vec3 &gvecJ_, const PS::F64vec3 &posI_, const PS::F64vec3 &posJ_, double phi0_)
        : gidI(gidI_), gidJ(gidJ_), globalIndexI(globalIndexI_), globalIndexJ(globalIndexJ_), phi0(phi0_), lambda(0) {
        // constructor
        gvecI = gvecI_;
        gvecJ = gvecJ_;
        posI = posI_;
        posJ = posJ_;
    }
};

class CalcNearForce {
  public:
    using colBlockThread_t = std::vector<std::deque<ColBlock> >;
    std::shared_ptr<colBlockThread_t> colBlockThread_ptr;

    // constructor
    CalcNearForce() {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        int totalThreads = omp_get_max_threads();
#else
        int totalThreads = 1;
#endif
        this->colBlockThread_ptr = std::make_shared<colBlockThread_t>();
        (*colBlockThread_ptr).resize(totalThreads);
        std::cout << "stress recoder size:" << (*colBlockThread_ptr).size() << std::endl;
    }

    // copy constructor
    CalcNearForce(const CalcNearForce &obj) : colBlockThread_ptr(obj.colBlockThread_ptr) {
        // copy the shared ptr
    }

    void operator()(const SphereNearForceEP *const ep_i, const PS::S32 Nip, const SphereNearForceEP *const ep_j,
                    const PS::S32 Njp, forceNearEP *const force) {

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        const int myThreadId = omp_get_thread_num();
#else
        const int myThreadId = 0;
#endif
        // int myRank;
        // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        for (PS::S32 i = 0; i < Nip; ++i) {
            force[i].clear();
            auto &sphereI = ep_i[i];
            double sepminI = 100000;
            PS::F64vec3 forceNearI = 0;
            PS::F64vec3 torqueNearI = 0;

            for (PS::S32 j = 0; j < Njp; ++j) {
                if (ep_j[j].gid == sphereI.gid) {
                    // no self force
                    continue;
                }
                auto &sphereJ = ep_j[j];
                auto rIJ = sphereJ.pos - sphereI.pos;
                const double rIJNorm = normVec3(rIJ);
                const double sep = rIJNorm - sphereI.radius - sphereJ.radius;
                sepminI = std::min(sepminI, sep);
                // process nearby force

                // save collision block
                // save only block gidI < gidJ, and gid
                if (sep < COLBUF * (sphereI.radius + sphereJ.radius) && sphereI.gid < sphereJ.gid) {
                    const PS::F64vec3 gvecI = rIJ * (-1.0 / rIJNorm);
                    const PS::F64vec3 gvecJ = -gvecI;
                    (*colBlockThread_ptr)[myThreadId].emplace_back(sphereI.gid, sphereJ.gid, sphereI.globalIndex,
                                                                   sphereJ.globalIndex, gvecI, gvecJ, sphereI.pos,
                                                                   sphereJ.pos, sep);
                }
            }
            force[i].forceNear = forceNearI;
            force[i].torqueNear = torqueNearI;
            force[i].sepmin = sepminI;
        }
    }
};

// spherical boundary.
class SphericalBoundary {
    double radius;
    PS::F64vec3 center;

  public:
    using colBlockThread_t = std::vector<std::deque<ColBlock> >;
    std::shared_ptr<colBlockThread_t> colBlockThread_ptr;
    // constructor
    SphericalBoundary() { initialize(0, PS::F64vec3(0, 0, 0)); }

    SphericalBoundary(double radius_, const PS::F64vec3 &center_) { initialize(radius_, center_); }

    void initialize(double radius_, const PS::F64vec3 &center_) {
        radius = radius_;
        center = center_;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        int totalThreads = omp_get_max_threads();
#else
        int totalThreads = 1;
#endif
        this->colBlockThread_ptr = std::make_shared<colBlockThread_t>();
        (*colBlockThread_ptr).resize(totalThreads);
        std::cout << "stress recoder size:" << (*colBlockThread_ptr).size() << std::endl;
        // only the first data slot of colBlock is used for sphere-boundary calculation
    }

    void clear() {
        for (auto &blockThread : *colBlockThread_ptr) {
            blockThread.clear(); // keep the number of threads but clear the contents
        }
    }

    void processParticle(const SphereFP &sphere) {
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
        const int threadId = omp_get_thread_num();
#else
        const int threadId = 0;
#endif
        processParticle(sphere, threadId);
    }

    void processParticle(const SphereFP &sphere, const int threadId) {
        const auto rvec = sphere.getPos() - center;
        const double rnorm = normVec3(rvec);
        const double sep = radius - (rnorm + sphere.radius);
        if (sep < sphere.radius * COLBUF) {
            auto gvec = -rvec * (1 / rnorm);
            (*colBlockThread_ptr)[threadId]
                .emplace_back(sphere.gid, -1, sphere.globalIndex, -1, gvec, -gvec, sphere.pos, center, sep);
            ;
        }
        return;
    }
};

class SphereSystem {

  public:
    const Config runConfig;
    int myRank;
    int nProcs;

    double sysTime;
    double snapTime;
    int id_snap;
    char dir_name[100];

    PS::ParticleSystem<SphereFP> systemSP;
    PS::DomainInfo dinfo;

    PS::TreeForForceShort<forceNearEP, SphereNearForceEP, SphereNearForceEP>::Gather treeNear;
    CalcNearForce calcNearForceFtr;
    /*
     * IMPORTANT: Scatter/Gather/Symmetry not only affects the RSearch.
     * It also changes the way that the EPI and EPJ are built.
     * eg. Here if treeFindBind is changed to Symmetry, the periodic images are not correctly found
     * */

    bool dumpxyz;
    bool dumpflow;

    std::unique_ptr<TRngPool> myRngPoolPtr;

    std::unique_ptr<HydroSphere::HydroSphereSystem> pointSphereSystemPtr;

    SphericalBoundary shell;

    bool withHydro; // an extra on\off of hydro

    void setInit(std::string);

    void moveAll(double);

    void calcSpringForce();
    void calcExtForce();

    void prepareT0();
    void stepForward();
    void output();
    void statistics();

    void overlapResolve(double);

    SphereSystem(std::string, std::string, int, char **); // constructor
    ~SphereSystem() = default;

    // delete copy constructor
    SphereSystem(const SphereSystem &) = delete;
    SphereSystem &operator=(const SphereSystem &) = delete;
};

//#pragma GCC pop_options

#endif
