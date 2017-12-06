#include "TpetraUtil.hpp"

void dumpTCMAT(const Teuchos::RCP<const TCMAT> &A, std::string filename) {
    filename = std::string("MatDump_") + filename;
#ifndef NDEBUG
    std::cout << "dumping" << filename << std::endl;
    // dump to matrixmarket format
    Tpetra::MatrixMarket::Writer<TCMAT> matDumper;
    matDumper.writeSparseFile(filename, A, filename, filename, true);
#endif
}

void dumpTMV(const Teuchos::RCP<const TMV> &A, std::string filename) {
    filename = std::string("VECDump_") + filename;
#ifndef NDEBUG
    std::cout << "dumping" << filename << std::endl;
    A->print(std::cout);
    // dump to matrixmarket format
    Tpetra::MatrixMarket::Writer<TMV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
#endif
}

void dumpTV(const Teuchos::RCP<const TV> &A, std::string filename) {
    filename = std::string("VECDump_") + filename;
#ifndef NDEBUG
    std::cout << "dumping" << filename << std::endl;
    A->print(std::cout);
    // dump to matrixmarket format
    Tpetra::MatrixMarket::Writer<TV> matDumper;
    matDumper.writeDenseFile(filename, A, filename, filename);
#endif
}
