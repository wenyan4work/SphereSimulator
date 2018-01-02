#include "TpetraUtil.hpp"

#ifdef DUMPTPETRA

void dumpTCMAT(const Teuchos::RCP<const TCMAT> &A, std::string filename) {
    filename = std::string("TCMAT_") + filename;
    std::cout << "dumping" << filename << std::endl;
    Tpetra::MatrixMarket::Writer<TCMAT> matDumper;
    matDumper.writeSparseFile(filename, A, filename, filename, true);
}

void dumpTMV(const Teuchos::RCP<const TMV> &A, std::string filename) {
    filename = std::string("TMV_") + filename;
    std::cout << "dumping" << filename << std::endl;
    A->print(std::cout);
    Tpetra::MatrixMarket::Writer<TMV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
}

void dumpTV(const Teuchos::RCP<const TV> &A, std::string filename) {
    filename = std::string("TV_") + filename;
    std::cout << "dumping" << filename << std::endl;
    A->print(std::cout);
    Tpetra::MatrixMarket::Writer<TV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
}

#endif
