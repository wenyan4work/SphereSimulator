#include "../../Buffer.hpp"
#include "../../SyncManager.hpp"

#include "mpi.h"
#include <cstdio>
#include <vector>

class Particle {
  public:
    int gid;
    std::vector<int> data;

    void Pack(std::vector<char> &buffer) const {
        Buffer mybuf(buffer);
        mybuf.pack(gid);
        mybuf.pack(data);
    }

    void Unpack(const std::vector<char> &buffer) {
        Buffer mybuf;
        mybuf.unpack(gid, buffer);
        mybuf.unpack(data, buffer);
    }

    void Reduce(Particle &other) {
        assert(gid == other.gid);
        std::copy(other.data.cbegin(), other.data.cend(), std::back_inserter(data));
    }

    void dump() {
        printf("--------------------\n");
        printf("%d,", gid);
        for (const auto &v : data) {
            printf("%d,", v);
        }
        printf("\n--------------------\n");
    }
};

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    int nProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // test pack and unpack
    if (rank == 0) {
        std::vector<char> buf;
        {
            Particle par;
            par.gid = 10;
            par.Pack(buf);
        }
        {
            Particle par2;
            par2.Unpack(buf);
            par2.dump();
        }
    }

    {

        std::vector<Particle> par;
        if (rank == 0) {
            par.resize(5);
            // fill data
            for (int i = 0; i < 5; i++) {
                par[i].gid = i;
                par[i].data.resize(i + 1);
            }
        }

        SyncManager<Particle> sync(&par);

        sync.sync();

        // display on rank 1
        if (rank == 1) {
            for (auto &obj : par) {
                obj.dump();
            }
        }

        // modify on every rank
        for (int p = 0; p < par.size(); p++) {
            auto &obj = par[p];
            obj.data.clear();
            for (int i = 0; i < p + 1; i++) {
                obj.data.push_back(rank);
            }
        }

        sync.reduce();
        printf("rank %d reduced\n", rank);

        // display on rank 0
        if (rank == 0) {
            for (auto &obj : par) {
                obj.dump();
            }
        }
    }

    MPI_Finalize();
    return 0;
}