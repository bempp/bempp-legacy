// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_discrete_sparse_boundary_operator_hpp
#define bempp_discrete_sparse_boundary_operator_hpp

#include "../common/common.hpp"
#include "bempp/common/config_ahmed.hpp"

#include "discrete_boundary_operator.hpp"
#include "../common/eigen_support.hpp"

#include "ahmed_aux_fwd.hpp"
#include "symmetry.hpp"
#include "transposition_mode.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp {

/** \ingroup discrete_boundary_operators
 *  \brief Discrete boundary operator stored as a sparse matrix.
 */
template <typename ValueType>
class DiscreteSparseBoundaryOperator
    : public DiscreteBoundaryOperator<ValueType> {
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

public:
  /** \brief Constructor.
   *
   *  \param[in] mat
   *    Sparse matrix that will be represented by the newly
   *    constructed operator. Must not be null.
   *  \param[in] symmetry
   *    Symmetry of the matrix. May be any combination of flags defined
   *    in the Symmetry enumeration type.
   *  \param[in] trans
   *    If different from NO_TRANSPOSE, the discrete operator will represent
   *    a transposed and/or complex-conjugated matrix \p mat. */
  DiscreteSparseBoundaryOperator(
      const shared_ptr<const RealSparseMatrix> &mat, int symmetry = NO_SYMMETRY,
      TranspositionMode trans = NO_TRANSPOSE);
private:
  DiscreteSparseBoundaryOperator();

public:

  virtual void dump() const;

  virtual Matrix<ValueType> asMatrix() const;

  virtual unsigned int rowCount() const;
  virtual unsigned int columnCount() const;

  virtual void addBlock(const std::vector<int> &rows,
                        const std::vector<int> &cols, const ValueType alpha,
                        Matrix<ValueType> &block) const;

  /** \brief Downcast a shared pointer to a DiscreteBoundaryOperator object to
   *  a shared pointer to a DiscreteSparseBoundaryOperator.
   *
   *  If the object referenced by \p discreteOperator is not in fact a
   *  DiscreteSparseBoundaryOperator, a std::bad_cast exception is thrown. */
  static shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>>
  castToSparse(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &
                   discreteOperator);


  /** \brief Return a shared pointer to the sparse matrix stored within
   *  this operator.
   *
   *  \note The discrete operator represents the matrix returned by this
   *  function *and possibly transposed and/or complex-conjugated*, depending on
   *  the value returned by transpositionMode(). */
  shared_ptr<const RealSparseMatrix> sparseMatrix() const;

  /** \brief Return the active sparse matrix transformation.
   *
   *  Indicates whether this operator represents the unmodified sparse matrix
   *  passed in the constructor or its transformation (transposition and/or
   *  conjugation). */
  TranspositionMode transpositionMode() const;

  /** \brief Return the symmetry type of the sparse matrix */
  inline int symmetryMode() const { return m_symmetry; }


private:
  /** \cond PRIVATE */
  virtual void applyBuiltInImpl(const TranspositionMode trans,
                                const Eigen::Ref<Vector<ValueType>> &x_in,
                                Eigen::Ref<Vector<ValueType>> y_inout,
                                const ValueType alpha,
                                const ValueType beta) const;
  bool isTransposed() const;

  /** \endcond */

private:
/** \cond PRIVATE */
  shared_ptr<const RealSparseMatrix> m_mat;
  int m_symmetry;
  TranspositionMode m_trans;
  /** \endcond */
};

} // namespace Bempp

#endif
