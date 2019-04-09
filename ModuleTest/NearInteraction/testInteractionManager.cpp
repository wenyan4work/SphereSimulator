// TODO:
// Test 1 species, partition, multiple interaction & multiple essential types
// Test 2 species, partition, multiple interaction & multiple essential types

#include "MPI/InteractionManager.hpp"
#include "Util/Buffer.hpp"

#include <random>

constexpr int DIM = 3;

constexpr int NPAR = 3000;
constexpr double radiusFrac = 0.05;

class FullParticle {
  public:
    int gid;
    double radius;
    double mass;
    double charge;
    double coord[DIM];
    double vel[DIM];
    double velGrav[DIM];
    double velElec[DIM];

    // required interface

    const double *Coord() const { return coord; }

    double Rad() const { return radius; }

    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        mybuff.pack(gid);
        mybuff.pack(radius);
        mybuff.pack(mass);
        mybuff.pack(charge);
        for (int i = 0; i < DIM; i++) {
            mybuff.pack(coord[i]);
            mybuff.pack(vel[i]);
            mybuff.pack(velGrav[i]);
            mybuff.pack(velElec[i]);
        }
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        mybuff.unpack(gid, buff);
        mybuff.unpack(radius, buff);
        mybuff.unpack(mass, buff);
        mybuff.unpack(charge, buff);
        for (int i = 0; i < DIM; i++) {
            mybuff.unpack(coord[i], buff);
            mybuff.unpack(vel[i], buff);
            mybuff.unpack(velGrav[i], buff);
            mybuff.unpack(velElec[i], buff);
        }
    }
};

class EssGrav {
  public:
    int gid;
    double mass;
    double radius;
    double coord[DIM];
    double velGrav[DIM] = {0, 0, 0};
    // required interface
    const double *Coord() const { return coord; }

    double Rad() const { return radius; }
    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        mybuff.pack(gid);
        mybuff.pack(radius);
        mybuff.pack(mass);
        for (int i = 0; i < DIM; i++) {
            mybuff.pack(coord[i]);
            mybuff.pack(velGrav[i]);
        }
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        mybuff.unpack(gid, buff);
        mybuff.unpack(radius, buff);
        mybuff.unpack(mass, buff);
        for (int i = 0; i < DIM; i++) {
            mybuff.unpack(coord[i], buff);
            mybuff.unpack(velGrav[i], buff);
        }
    }

    void CopyFromFull(const FullParticle &fp) {
        gid = fp.gid;
        mass = fp.mass;
        radius = fp.radius; // ess can have different radius from fp. here it is simplified
        for (int i = 0; i < DIM; i++) {
            coord[i] = fp.coord[i];
        }
    }

    void dump() const { std::cout << velGrav[0] << " " << velGrav[1] << " " << velGrav[2] << std::endl; }
};

class EssElec {
  public:
    int gid;
    double charge;
    double radius;
    double coord[DIM];
    double velElec[DIM] = {0, 0, 0};
    // required interface
    inline const double *Coord() const { return coord; }

    inline double Rad() const { return radius; }

    void Pack(std::vector<char> &buff) const {
        Buffer mybuff(buff);
        mybuff.pack(gid);
        mybuff.pack(radius);
        mybuff.pack(charge);
        for (int i = 0; i < DIM; i++) {
            mybuff.pack(coord[i]);
            mybuff.pack(velElec[i]);
        }
    }

    void Unpack(const std::vector<char> &buff) {
        Buffer mybuff;
        mybuff.unpack(gid, buff);
        mybuff.unpack(radius, buff);
        mybuff.unpack(charge, buff);
        for (int i = 0; i < DIM; i++) {
            mybuff.unpack(coord[i], buff);
            mybuff.unpack(velElec[i], buff);
        }
    }

    void CopyFromFull(const FullParticle &fp) {
        gid = fp.gid;
        charge = fp.charge;
        radius = fp.radius * 0.3; // WARNING:try a different radius to test
        for (int i = 0; i < DIM; i++) {
            coord[i] = fp.coord[i];
        }
    }

    void dump() const { std::cout << velElec[0] << " " << velElec[1] << " " << velElec[2] << std::endl; }
};

class Gravity {
  public:
    inline void operator()(EssGrav &t, const EssGrav &s, const std::array<double, DIM> &srcShift) {
#ifdef DEBUGINTERACT
        std::cout << t.gid << " " << s.gid << std::endl;
        std::cout << t.coord[0] << " " << s.coord[0] << " " << srcShift[0] << std::endl;
        std::cout << t.coord[1] << " " << s.coord[1] << " " << srcShift[1] << std::endl;
        std::cout << t.coord[2] << " " << s.coord[2] << " " << srcShift[2] << std::endl;
#endif
        double dr[DIM];
        double rnorm = 0;
        for (int i = 0; i < DIM; i++) {
            dr[i] = t.coord[i] - (s.coord[i] + srcShift[i]);
            rnorm += dr[i] * dr[i];
        }
        rnorm += 0.001; // avoid div0 error
        rnorm = sqrt(rnorm);
        for (int i = 0; i < DIM; i++) {
            t.velGrav[i] += s.mass * t.mass * dr[i] / pow(rnorm, 3);
        }
    }
};

class Elec {
  public:
    inline void operator()(EssElec &t, const EssElec &s, const std::array<double, DIM> &srcShift) {
#ifdef DEBUGINTERACT
        std::cout << t.gid << " " << s.gid << std::endl;
        std::cout << t.coord[0] << " " << s.coord[0] << " " << srcShift[0] << std::endl;
        std::cout << t.coord[1] << " " << s.coord[1] << " " << srcShift[1] << std::endl;
        std::cout << t.coord[2] << " " << s.coord[2] << " " << srcShift[2] << std::endl;
#endif
        double dr[DIM];
        double rnorm = 0;
        for (int i = 0; i < DIM; i++) {
            dr[i] = t.coord[i] - (s.coord[i] + srcShift[i]);
            rnorm += dr[i] * dr[i];
        }
        rnorm += 0.001; // avoid div0 error
        rnorm = sqrt(rnorm);
        for (int i = 0; i < DIM; i++) {
            t.velElec[i] += s.charge * t.charge * dr[i] / pow(rnorm, 3);
        }
    }
};

void testOneSpecies(const int NPAR) {
    int myRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // initialize to random
    std::vector<FullParticle> particle;

    // std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> dis(0, 1);

    std::array<double, 3> boxLow = {10.0, 10.0, 10.0};
    std::array<double, 3> box = {10.0, 20.0, 30.0};
    std::array<double, 3> boxHigh;

    boxHigh[0] = boxLow[0] + box[0];
    boxHigh[1] = boxLow[1] + box[1];
    boxHigh[2] = boxLow[2] + box[2];

    if (myRank == 0) {
        particle.resize(NPAR);
        for (int i = 0; i < NPAR; i++) {
            auto &p = particle[i];
            p.gid = myRank * NPAR + i;
            for (int j = 0; j < DIM; j++) {
                p.coord[j] = dis(gen) * box[j];
            }
            p.mass = 1;
            p.charge = 1;
            p.radius = radiusFrac * box[0];
            for (int j = 0; j < DIM; j++) {
                p.vel[j] = 0;
                p.velGrav[j] = 0;
                p.velElec[j] = 0;
            }
        }
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2]);
    // }

    // initialize
    InteractionManager<double, DIM, FullParticle, FullParticle> interactManager(&particle, &particle);

    // create a nearInteractor for full particle
    auto nearInteractFullParPtr = interactManager.getNewNearInteraction();
    interactManager.partitionObject(nearInteractFullParPtr);
    printf("partition complete\n");

    auto nearInteractPtr = interactManager.getNewNearInteraction();

    // near interactor for Grav
    std::vector<EssGrav> srcEss, trgEss;
    Gravity grav;
    interactManager.setupEssVec(srcEss, trgEss);
    fprintf(stderr, "setupEssVec complete\n");
    interactManager.setupNearInteractor(nearInteractPtr, srcEss, trgEss);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssGrav, EssGrav, Gravity>(nearInteractPtr, srcEss, trgEss, grav);
    fprintf(stderr, "calcNearInteraction complete\n");

    // write back
    assert(trgEss.size() == particle.size());
    const int npar = particle.size();
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velGrav[0] = trgEss[i].velGrav[0];
        particle[i].velGrav[1] = trgEss[i].velGrav[1];
        particle[i].velGrav[2] = trgEss[i].velGrav[2];
    }

    nearInteractPtr->Barrier();

    // near interactor for Elec
    std::vector<EssElec> srcEssElec, trgEssElec;
    Elec elec;
    interactManager.setupEssVec(srcEssElec, trgEssElec);
    fprintf(stderr, "setupEssVec complete\n");
    nearInteractPtr = interactManager.getNewNearInteraction();
    interactManager.setupNearInteractor(nearInteractPtr, srcEssElec, trgEssElec);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssElec, EssElec, Elec>(nearInteractPtr, srcEssElec, trgEssElec, elec);
    fprintf(stderr, "calcNearInteraction complete\n");
    assert(trgEssElec.size() == particle.size());
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velElec[0] = trgEssElec[i].velElec[0];
        particle[i].velElec[1] = trgEssElec[i].velElec[1];
        particle[i].velElec[2] = trgEssElec[i].velElec[2];
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f, %f,%f,%f,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2], p.velGrav[0],
    //            p.velGrav[1], p.velGrav[2], p.velElec[0], p.velElec[1], p.velElec[2]);
    // }

    // check sum zero
    double vex = 0, vey = 0, vez = 0;
    double vgx = 0, vgy = 0, vgz = 0;
#pragma omp parallel for reduction(+ : vex, vey, vez)
    for (int i = 0; i < npar; i++) {
        vex += particle[i].velElec[0];
        vey += particle[i].velElec[1];
        vez += particle[i].velElec[2];
    }
    printf("sum velElec: %f,%f,%f\n", vex, vey, vez);

#pragma omp parallel for reduction(+ : vgx, vgy, vgz)
    for (int i = 0; i < npar; i++) {
        vgx += particle[i].velGrav[0];
        vgy += particle[i].velGrav[1];
        vgz += particle[i].velGrav[2];
    }
    printf("sum velGrav: %f,%f,%f\n", vgx, vgy, vgz);

    {
        double vGlobalSum[3] = {0, 0, 0};
        double ve[3] = {vex, vey, vez};
        MPI_Reduce(ve, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velElec: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    {
        double vGlobalSum[3] = {0, 0, 0};
        double vg[3] = {vgx, vgy, vgz};
        MPI_Reduce(vg, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velGrav: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    return;
}

void testOneSpeciesPBCX(const int NPAR) {
    int myRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // initialize to random
    std::vector<FullParticle> particle;

    // std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> dis(0, 1);

    std::array<bool, 3> pbcFLag = {1, 0, 0};
    std::array<double, 3> boxLow = {10.0, 10.0, 10.0};
    std::array<double, 3> box = {10.0, 20.0, 30.0};
    std::array<double, 3> boxHigh;

    boxHigh[0] = boxLow[0] + box[0];
    boxHigh[1] = boxLow[1] + box[1];
    boxHigh[2] = boxLow[2] + box[2];

    if (myRank == 0) {
        particle.resize(NPAR);
        for (int i = 0; i < NPAR; i++) {
            auto &p = particle[i];
            p.gid = myRank * NPAR + i;
            for (int j = 0; j < DIM; j++) {
                p.coord[j] = dis(gen) * box[j];
            }
            p.mass = 1;
            p.charge = 1;
            p.radius = radiusFrac * box[0];
            for (int j = 0; j < DIM; j++) {
                p.vel[j] = 0;
                p.velGrav[j] = 0;
                p.velElec[j] = 0;
            }
        }
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2]);
    // }

    // initialize
    InteractionManager<double, DIM, FullParticle, FullParticle> interactManager(&particle, &particle);
    interactManager.setPBCBox(pbcFLag, boxLow, boxHigh);

    // create a nearInteractor for full particle
    auto nearInteractFullParPtr = interactManager.getNewNearInteraction();
    interactManager.partitionObject(nearInteractFullParPtr);
    printf("partition complete\n");

    auto nearInteractPtr = interactManager.getNewNearInteraction();

    // near interactor for Grav
    std::vector<EssGrav> srcEss, trgEss;
    Gravity grav;
    interactManager.setupEssVec(srcEss, trgEss);
    fprintf(stderr, "setupEssVec complete\n");
    interactManager.setupNearInteractor(nearInteractPtr, srcEss, trgEss);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssGrav, EssGrav, Gravity>(nearInteractPtr, srcEss, trgEss, grav);
    fprintf(stderr, "calcNearInteraction complete\n");

    // write back
    assert(trgEss.size() == particle.size());
    const int npar = particle.size();
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velGrav[0] = trgEss[i].velGrav[0];
        particle[i].velGrav[1] = trgEss[i].velGrav[1];
        particle[i].velGrav[2] = trgEss[i].velGrav[2];
    }

    nearInteractPtr->Barrier();

    // near interactor for Elec
    std::vector<EssElec> srcEssElec, trgEssElec;
    Elec elec;
    interactManager.setupEssVec(srcEssElec, trgEssElec);
    fprintf(stderr, "setupEssVec complete\n");
    nearInteractPtr = interactManager.getNewNearInteraction();
    interactManager.setupNearInteractor(nearInteractPtr, srcEssElec, trgEssElec);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssElec, EssElec, Elec>(nearInteractPtr, srcEssElec, trgEssElec, elec);
    fprintf(stderr, "calcNearInteraction complete\n");
    assert(trgEssElec.size() == particle.size());
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velElec[0] = trgEssElec[i].velElec[0];
        particle[i].velElec[1] = trgEssElec[i].velElec[1];
        particle[i].velElec[2] = trgEssElec[i].velElec[2];
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f, %f,%f,%f,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2], p.velGrav[0],
    //            p.velGrav[1], p.velGrav[2], p.velElec[0], p.velElec[1], p.velElec[2]);
    // }

    // check sum
    double vex = 0, vey = 0, vez = 0;
    double vgx = 0, vgy = 0, vgz = 0;
#pragma omp parallel for reduction(+ : vex, vey, vez)
    for (int i = 0; i < npar; i++) {
        vex += particle[i].velElec[0];
        vey += particle[i].velElec[1];
        vez += particle[i].velElec[2];
    }
    printf("sum velElec: %f,%f,%f\n", vex, vey, vez);

#pragma omp parallel for reduction(+ : vgx, vgy, vgz)
    for (int i = 0; i < npar; i++) {
        vgx += particle[i].velGrav[0];
        vgy += particle[i].velGrav[1];
        vgz += particle[i].velGrav[2];
    }
    printf("sum velGrav: %f,%f,%f\n", vgx, vgy, vgz);

    {
        double vGlobalSum[3] = {0, 0, 0};
        double ve[3] = {vex, vey, vez};
        MPI_Reduce(ve, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velElec: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    {
        double vGlobalSum[3] = {0, 0, 0};
        double vg[3] = {vgx, vgy, vgz};
        MPI_Reduce(vg, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velGrav: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    return;
}

void testOneSpeciesPBCYZ(const int NPAR) {
    int myRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // initialize to random
    std::vector<FullParticle> particle;

    // std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> dis(0, 1);

    std::array<bool, 3> pbcFLag = {0, 1, 1};
    std::array<double, 3> boxLow = {10.0, 10.0, 10.0};
    std::array<double, 3> box = {10.0, 20.0, 30.0};
    std::array<double, 3> boxHigh;

    boxHigh[0] = boxLow[0] + box[0];
    boxHigh[1] = boxLow[1] + box[1];
    boxHigh[2] = boxLow[2] + box[2];

    if (myRank == 0) {
        particle.resize(NPAR);
        for (int i = 0; i < NPAR; i++) {
            auto &p = particle[i];
            p.gid = myRank * NPAR + i;
            for (int j = 0; j < DIM; j++) {
                p.coord[j] = dis(gen) * box[j];
            }
            p.mass = 1;
            p.charge = 1;
            p.radius = radiusFrac * box[0];
            for (int j = 0; j < DIM; j++) {
                p.vel[j] = 0;
                p.velGrav[j] = 0;
                p.velElec[j] = 0;
            }
        }
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2]);
    // }

    // initialize
    InteractionManager<double, DIM, FullParticle, FullParticle> interactManager(&particle, &particle);
    interactManager.setPBCBox(pbcFLag, boxLow, boxHigh);

    // create a nearInteractor for full particle
    auto nearInteractFullParPtr = interactManager.getNewNearInteraction();
    interactManager.partitionObject(nearInteractFullParPtr);
    printf("partition complete\n");

    auto nearInteractPtr = interactManager.getNewNearInteraction();

    // near interactor for Grav
    std::vector<EssGrav> srcEss, trgEss;
    Gravity grav;
    interactManager.setupEssVec(srcEss, trgEss);
    fprintf(stderr, "setupEssVec complete\n");
    interactManager.setupNearInteractor(nearInteractPtr, srcEss, trgEss);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssGrav, EssGrav, Gravity>(nearInteractPtr, srcEss, trgEss, grav);
    fprintf(stderr, "calcNearInteraction complete\n");

    // write back
    assert(trgEss.size() == particle.size());
    const int npar = particle.size();
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velGrav[0] = trgEss[i].velGrav[0];
        particle[i].velGrav[1] = trgEss[i].velGrav[1];
        particle[i].velGrav[2] = trgEss[i].velGrav[2];
    }

    nearInteractPtr->Barrier();

    // near interactor for Elec
    std::vector<EssElec> srcEssElec, trgEssElec;
    Elec elec;
    interactManager.setupEssVec(srcEssElec, trgEssElec);
    fprintf(stderr, "setupEssVec complete\n");
    nearInteractPtr = interactManager.getNewNearInteraction();
    interactManager.setupNearInteractor(nearInteractPtr, srcEssElec, trgEssElec);
    fprintf(stderr, "setupNearInteractor complete\n");
    interactManager.calcNearInteraction<EssElec, EssElec, Elec>(nearInteractPtr, srcEssElec, trgEssElec, elec);
    fprintf(stderr, "calcNearInteraction complete\n");
    assert(trgEssElec.size() == particle.size());
#pragma omp parallel for
    for (int i = 0; i < npar; i++) {
        particle[i].velElec[0] = trgEssElec[i].velElec[0];
        particle[i].velElec[1] = trgEssElec[i].velElec[1];
        particle[i].velElec[2] = trgEssElec[i].velElec[2];
    }

    // // output
    // for (const auto &p : particle) {
    //     printf("%d,%f,%f,%f, %f,%f,%f,%f,%f,%f\n", p.gid, p.coord[0], p.coord[1], p.coord[2], p.velGrav[0],
    //            p.velGrav[1], p.velGrav[2], p.velElec[0], p.velElec[1], p.velElec[2]);
    // }

    // check sum
    double vex = 0, vey = 0, vez = 0;
    double vgx = 0, vgy = 0, vgz = 0;
#pragma omp parallel for reduction(+ : vex, vey, vez)
    for (int i = 0; i < npar; i++) {
        vex += particle[i].velElec[0];
        vey += particle[i].velElec[1];
        vez += particle[i].velElec[2];
    }
    printf("sum velElec: %f,%f,%f\n", vex, vey, vez);

#pragma omp parallel for reduction(+ : vgx, vgy, vgz)
    for (int i = 0; i < npar; i++) {
        vgx += particle[i].velGrav[0];
        vgy += particle[i].velGrav[1];
        vgz += particle[i].velGrav[2];
    }
    printf("sum velGrav: %f,%f,%f\n", vgx, vgy, vgz);

    {
        double vGlobalSum[3] = {0, 0, 0};
        double ve[3] = {vex, vey, vez};
        MPI_Reduce(ve, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velElec: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    {
        double vGlobalSum[3] = {0, 0, 0};
        double vg[3] = {vgx, vgy, vgz};
        MPI_Reduce(vg, vGlobalSum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            printf("Global sum velGrav: %g,%g,%g\n", vGlobalSum[0], vGlobalSum[1], vGlobalSum[2]);
    }

    return;
}

int main(int argc, char **argv) {
    // omp_set_num_threads(1);

    MPI_Init(&argc, &argv);

    printf("testing free space BC\n");
    testOneSpecies(NPAR);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("testing PBC X\n");
    testOneSpeciesPBCX(NPAR);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("testing PBC YZ\n");
    testOneSpeciesPBCYZ(NPAR);

    MPI_Finalize();

    return 0;
}
