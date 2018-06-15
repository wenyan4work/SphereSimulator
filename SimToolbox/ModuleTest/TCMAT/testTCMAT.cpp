#include "../../Trilinos/TpetraUtil.hpp"
#include <mpi.h>

#include <vector>

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    Teuchos::RCP<const TCOMM> commRcp = getMPIWORLDTCOMM();

    std::vector<int> gColIndexOnLocal = {1 + commRcp->getRank(), 3 + commRcp->getRank()};
    std::vector<int> gRowIndexOnLocal = {commRcp->getRank() + 1};

    Teuchos::RCP<TMAP> rowMapRcp = getTMAPFromLocalSize(1, commRcp);
    Teuchos::RCP<TMAP> colOpMapRcp = getTMAPFromLocalSize(gColIndexOnLocal.size(), commRcp);
    Teuchos::RCP<TMAP> colMapRcp = Teuchos::rcp<TMAP>(new TMAP(5, gColIndexOnLocal.data(), 2, 0, commRcp));

    // entries
    Kokkos::View<size_t *> rowPointers("rowPointers", gRowIndexOnLocal.size() + 1);
    rowPointers[0] = 0;
    rowPointers[1] = rowPointers[0] + gColIndexOnLocal.size();
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[1]);
    Kokkos::View<double *> values("values", rowPointers[1]);
    columnIndices[0] = gColIndexOnLocal[0];
    columnIndices[1] = gColIndexOnLocal[1];
    values[0] = gColIndexOnLocal[0];
    values[1] = gColIndexOnLocal[1];

    auto &colmap = *colMapRcp;
    const int colIndexCount = gColIndexOnLocal.size();
#pragma omp parallel for
    for (int i = 0; i < colIndexCount; i++) {
        columnIndices[i] = colmap.getLocalElement(columnIndices[i]);
    }

    Teuchos::RCP<TCMAT> matTestRcp = Teuchos::rcp(new TCMAT(rowMapRcp, colMapRcp, rowPointers, columnIndices, values));
    matTestRcp->fillComplete(colOpMapRcp, rowMapRcp); // domainMap, rangeMap

    dumpTCMAT(matTestRcp, "matTestRcp");

    MPI_Finalize();
    return 0;
}