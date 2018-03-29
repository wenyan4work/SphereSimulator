/*
 * TpetraDef.hpp
 *  Created on: Dec 20, 2016
 *      Author: wyan
 */

#ifndef TPETRAUTIL_HPP_
#define TPETRAUTIL_HPP_

// Utility
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_oblackholestream.hpp>

// Container
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

// Solver
#include <BelosOperatorTraits.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

// Preconditioner
#include <Ifpack2_Factory.hpp>

// no need to specify node type for new version of Tpetra. It defaults to
// Kokkos::default, which is openmp
// typedef Tpetra::Details::DefaultTypes::node_type TNODE;
using TCOMM = Teuchos::Comm<int>;
using TMAP = Tpetra::Map<int, int>;
using TOP = Tpetra::Operator<double, int, int>;
using TCMAT = Tpetra::CrsMatrix<double, int, int>;
using TMV = Tpetra::MultiVector<double, int, int>;
using TV = Tpetra::Vector<double, int, int>;

// explicit full specialization of OperatorTraits for Tpetra objects.
namespace Belos {
template <>
class OperatorTraits<::TOP::scalar_type, ::TMV, ::TOP> {
  public:
    static void Apply(const ::TOP &Op, const ::TMV &X, ::TMV &Y, const ETrans trans = NOTRANS) {
        Teuchos::ETransp teuchosTrans = Teuchos::NO_TRANS;
        if (trans == NOTRANS) {
            teuchosTrans = Teuchos::NO_TRANS;
        } else if (trans == TRANS) {
            teuchosTrans = Teuchos::TRANS;
        } else if (trans == CONJTRANS) {
            teuchosTrans = Teuchos::CONJ_TRANS;
        } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                       "Belos::OperatorTraits::Apply: Invalid "
                                       "'trans' value "
                                           << trans << ".  Valid values are NOTRANS=" << NOTRANS << ", TRANS=" << TRANS
                                           << ", and CONJTRANS=" << CONJTRANS << ".");
        }
        Op.apply(X, Y, teuchosTrans);
    }

    static bool HasApplyTranspose(const ::TOP &Op) { return Op.hasTransposeApply(); }
};
} // namespace Belos

// utility functions
void dumpTCMAT(const Teuchos::RCP<const TCMAT> &A, std::string filename);
void dumpTMV(const Teuchos::RCP<const TMV> &A, std::string filename);
void dumpTV(const Teuchos::RCP<const TV> &A, std::string filename);

// the default TCOMM corresponding to MPI_COMM_WORLD
Teuchos::RCP<const TCOMM> getMPIWORLDTCOMM();

// return a fully copied TMAP with a given global size
Teuchos::RCP<TMAP> getFullCopyTMAPFromGlobalSize(const int &globalSize, Teuchos::RCP<const TCOMM> &commRcp);
// return a contiguous TMAP from local Size
Teuchos::RCP<TMAP> getTMAPFromLocalSize(const int &localSize, Teuchos::RCP<const TCOMM> &commRcp);

// contiguous TV init from a vector
Teuchos::RCP<TV> getTVFromVector(const std::vector<double> &in, Teuchos::RCP<const TCOMM> &commRcp);

// contiguous TMV init from a vector of vector. localsize= min_k in[k].size()
Teuchos::RCP<TMV> getTMVFromVector(const std::vector<std::vector<double>> &in, Teuchos::RCP<const TCOMM> &commRcp);


#endif /* TPETRAUTIL_HPP_ */
