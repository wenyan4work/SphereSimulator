/*
 * TpetraDef.hpp
 *
 *  Created on: Dec 20, 2016
 *      Author: wyan
 */

#ifndef TPETRAUTIL_HPP_
#define TPETRAUTIL_HPP_

#include <BelosOperatorTraits.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Ifpack2_Factory.hpp>

#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

#include <MatrixMarket_Tpetra.hpp>

// no need to specify node type for new version of Tpetra. It defaults to
// Kokkos::default, which is openmp
// typedef Tpetra::Details::DefaultTypes::node_type TNODE;
typedef Tpetra::Map<int, int> TMAP;
typedef Tpetra::CrsMatrix<double, int, int> TCMAT;
typedef Tpetra::Operator<double, int, int> TOP;
typedef Tpetra::MultiVector<double, int, int> TMV;
typedef Tpetra::Vector<double, int, int> TV;

typedef Teuchos::Comm<int> TCOMM;

void dumpTCMAT(const Teuchos::RCP<const TCMAT> &A, std::string filename);

void dumpTMV(const Teuchos::RCP<const TMV> &A, std::string filename);

void dumpTV(const Teuchos::RCP<const TV> &A, std::string filename);

// explicit full specialization of OperatorTraits for Tpetra objects.
namespace Belos {
    template<>
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
                        "Belos::OperatorTraits::Apply: Invalid " "'trans' value " << trans << ".  Valid values are NOTRANS=" << NOTRANS << ", TRANS=" << TRANS << ", and CONJTRANS=" << CONJTRANS << ".");
            }
            Op.apply(X, Y, teuchosTrans);
        }

        static bool HasApplyTranspose(const ::TOP &Op) {
            return Op.hasTransposeApply();
        }
    };
}

#endif /* TPETRAUTIL_HPP_ */
