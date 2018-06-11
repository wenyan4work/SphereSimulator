
#ifndef ZDD_HPP_
#define ZDD_HPP_

#include <vector>

#include "mpi.h"
#include <zoltan_dd_cpp.h>
#include <zoltan_types.h>

template <class DATA_TYPE> class ZDD {
  public:
    typedef ZOLTAN_ID_TYPE ID_TYPE;

    int rankSize;
    int myRank;

    Zoltan_DD findZDD;

    std::vector<ID_TYPE> findID;
    std::vector<DATA_TYPE> findData;
    std::vector<ID_TYPE> localID;
    std::vector<DATA_TYPE> localData;

    // constructor, with an estimate of buffer list size
    ZDD(int nEst) {
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        MPI_Comm_size(MPI_COMM_WORLD, &rankSize);

        this->findZDD.Create(MPI_COMM_WORLD, 1, 0, sizeof(DATA_TYPE), 0, 0);
#ifdef ZDDDEBUG
        this->findZDD.Create(MPI_COMM_WORLD, 1, 0, sizeof(DATA_TYPE), 0, 9); // debug verbose level 9
        findZDD.Stats();
#endif
        findID.reserve(nEst);
        findData.reserve(nEst);
        localID.reserve(nEst);
        localData.reserve(nEst);
    }

    // manual clear up
    void clearAll() {
        findID.clear();
        findData.clear();
        localID.clear();
        localData.clear();
    }

    //  destructor
    ~ZDD() {
        MPI_Barrier(MPI_COMM_WORLD);
//	this->nbFindZDD.~Zoltan_DD();
#ifdef ZDDDEBUG
        std::cout << "ZDD destructed" << std::endl;
#endif
    }

    // must have an effective nblocal data list
    int buildIndex() {

#ifdef ZDDDEBUG
        printf("zoltan print 1 buildIndex\n");
        this->findZDD.Print();
#endif

        auto *idPtr = this->localID.data();
        auto *dataPtr = this->localData.data();
        int error;
        //	error = nbFindZDD.Remove(idPtr, nbLocalID.size()); // remove the old info
        error = findZDD.Update(idPtr, NULL, (char *)dataPtr, NULL, localID.size());

#ifdef ZDDDEBUG
        printf("zoltan print 2 fildIndex\n");
        findZDD.Print();
        printf("%d\n", error);
        findZDD.Stats();
#endif

        return error;
    }

    // must have valid findID list
    int find() {
        if (findData.size() < findID.size()) {
            findData.resize(findID.size()); // make sure size match
        }
        auto *idPtr = findID.data();
        auto *dataPtr = findData.data();
        // ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid, char *data, int *part, const int &
        // count, int *owner
        int status = findZDD.Find(idPtr, NULL, (char *)dataPtr, NULL, findID.size(), NULL);
        return status;
    }
};

#endif /* ZDD_HPP_ */
