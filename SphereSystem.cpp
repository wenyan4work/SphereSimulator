#include <cassert>
#include <cmath>
#include <vector>

#include <mpi.h>
#include <omp.h>

#include "SphereSystem.hpp"
#include "Util/IOHelper.hpp"

SphereSystem::SphereSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv)
    : runConfig(configFile) {
    int mpiflag;
    MPI_Initialized(&mpiflag);
    assert(mpiflag);

    commRcp = getMPIWORLDTCOMM();

    sysTime = 0;
    snapTime = 0;
    snapID = 0;

    // TRNG pool must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(runConfig.rngSeed);

    if (commRcp->getRank() == 0) {
        // prepare the output directory
        IOHelper::makeSubFolder("./result");
        // initialize on head rank
        setInitial(posFile);
    } else {
        sphere.reserve(runConfig.sphereNumber / commRcp->getSize() * 2);
    }

    // initial exchange of particles
    MPI_Barrier(MPI_COMM_WORLD);
    partition();
    MPI_Barrier(MPI_COMM_WORLD);

    printf("local: %d spheres on process %d\n", sphere.size(), commRcp->getRank());

    // initialize force calculate tree
    // initial resolve and reset move to zero
    // 5 pass to resolve initial overlaps

    printf("SphereSystem initialized on process: %d \n", commRcp->getRank());
}

void SphereSystem::setInitial(const std::string &initPos) {}

bool SphereSystem::readXYZ(const std::string &filename) {
    std::ifstream myfile(filename);
    if (!myfile.is_open()) {
        printf("Error: cannot open initial position file. use random initial configuration\n");
        return false;
    }
    std::string line;
    getline(myfile, line); //
    getline(myfile, line); // read two header lines

    while (std::getline(myfile, line)) {
        // read 5 numbers: gid, radius, x,y,z
        // set radiusCollision = radius
        std::istringstream liness(line);
        int gid;
        char type;
        double px, py, pz;
        double radius;
        liness >> type >> gid >> radius >> px >> py >> pz;
        std::cout << gid << " " << radius << " " << px << " " << py << " "
                  << " " << pz << "\n";
        sphere.emplace_back(gid, radius, radius * runConfig.sphereRadiusCollisionRatio, Evec3(px, py, pz));
    }
    myfile.close();

    std::cout << "read initial finished" << std::endl;
    return true;
}

void SphereSystem::writeXYZ(const std::string &baseFolder) {
    // each rank write a xyz file separately. rank 0 write the two lines of header
    std::string name = baseFolder + std::string("SphereXYZ_") + std::to_string(snapID) + std::string("_r") +
                       std::to_string(commRcp->getRank()) + "_.xyz";
    FILE *fptr = fopen(name.c_str(), "w");
    const int nSphereLocal = sphere.size();
    int nSphere = nSphereLocal;
    if (commRcp->getSize() > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &nSphere, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    if (commRcp->getRank() == 0) {
        fprintf(fptr, "%d\n", nSphere);
        fprintf(fptr, "SphereNumber: %d, time: %lf\n", nSphere, sysTime);
    }

    for (const auto &s : sphere) {
        fprintf(fptr, "S\t%8d\t%.6g\t%.6g\t%.6g\t%.6g\n", s.gid, s.radius, s.pos[0], s.pos[1], s.pos[2]);
    }

    fclose(fptr);
}

void SphereSystem::writeVTK(const std::string &baseFolder) {
    Sphere::writeVTP(sphere, baseFolder, std::to_string(snapID), commRcp->getRank());
    Sphere::writeVTU(sphere, baseFolder, std::to_string(snapID), commRcp->getRank());
    // write parallel head
    if (commRcp->getSize() > 1 && commRcp->getRank() == 0) {
        Sphere::writePVTP(baseFolder, "000", commRcp->getSize());
        std::vector<std::pair<int, std::string>> dataFields;
        std::vector<IOHelper::IOTYPE> types = {IOHelper::IOTYPE::Float64, IOHelper::IOTYPE::Float64};
        for (auto &layer : sphere[0].sphLayer) {
            dataFields.emplace_back(std::pair<int, std::string>(layer.second->dimension, layer.second->name));
        }
        dataFields.emplace_back(std::pair<int, std::string>(3, "stk"));
        Sphere::writePVTU(dataFields, types, baseFolder, std::to_string(snapID), commRcp->getSize());
    }
}

void SphereSystem::output() {
    // every 200 snaps in a sub folder
    const int num = 200;
    int k = snapID / num;
    int low = k * num, high = k * num + num - 1;
    std::string baseFolder =
        "result/" + std::to_string(low) + std::string("-") + std::to_string(high) + std::string("/");
    IOHelper::makeSubFolder(baseFolder);
    writeXYZ(baseFolder);
    writeVTK(baseFolder);
}