#ifndef bempp_sparse_cholesky_hpp
#define bempp_sparse_cholesky_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

class Epetra_CrsMatrix;

namespace Bempp {

/** \relates DiscreteSparseBoundaryOperator
 *  \brief Return the first factor of the Cholesky factorisation of matrix \p
 *mat.
 *
 *  The matrix \p mat is assumed to be symmetric, positive-definite and
 *  block-diagonal (effectively it should be the mass matrix of a space of
 *  discontinuous functions). The routine returns a shared pointer to a matrix
 *  L such that
 *
 *  \f[ A = L L^T, \f]
 *
 *  where \f$A\f$ is the matrix represented by \p mat.
 */
shared_ptr<Epetra_CrsMatrix> sparseCholesky(const Epetra_CrsMatrix &mat);

} // namespace Bempp

#endif
