#ifndef PRECOND_HPP_
#define PRECOND_HPP_

#include "TpetraUtil.hpp"

Teuchos::RCP<TOP> createILUTPreconditioner(const Teuchos::RCP<const TCMAT> &, double, double);

Teuchos::RCP<TOP> createPlnPreconditioner(const Teuchos::RCP<const TCMAT> &);

Teuchos::RCP<TOP> createKinvPreconditioner(const Teuchos::RCP<const TCMAT> &, const Teuchos::RCP<const TMV> &);
#endif