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

#ifndef bempp_aca_approximate_lu_inverse_hpp
#define bempp_aca_approximate_lu_inverse_hpp

#include "../common/common.hpp"

#include "bempp/common/config_trilinos.hpp"
#include "discrete_boundary_operator.hpp"

#include "../fiber/scalar_traits.hpp"
#include "../fiber/verbosity_level.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

using Fiber::VerbosityLevel;

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteAcaBoundaryOperator;
/** \endcond */

/** \ingroup composite_discrete_operators
 *  \brief Approximate LU decomposition of a H-matrix
 */
template <typename ValueType>
class AcaApproximateLuInverse : public DiscreteBoundaryOperator<ValueType> {
public:
  typedef typename Fiber::ScalarTraits<ValueType>::RealType MagnitudeType;

  /** \brief Construct an approximate LU decomposition of a H-matrix.

  \param[in] fwdOp  Operator represented internally as a H-matrix.
  \param[in] delta  Requested approximation accuracy (M. Bebendorf recommends
                    delta = 0.1). */
  AcaApproximateLuInverse(const DiscreteAcaBoundaryOperator<ValueType> &fwdOp,
                          MagnitudeType delta,
                          VerbosityLevel::Level verbosityLevel =
                              VerbosityLevel::DEFAULT);

  virtual ~AcaApproximateLuInverse();

  virtual unsigned int rowCount() const;
  virtual unsigned int columnCount() const;

  virtual void addBlock(const std::vector<int> &rows,
                        const std::vector<int> &cols, const ValueType alpha,
                        Matrix<ValueType> &block) const;

#ifdef WITH_TRILINOS
public:
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>> domain() const;
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>> range() const;

protected:
  virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
#endif

private:
  virtual void applyBuiltInImpl(const TranspositionMode trans,
                                const Vector<ValueType> &x_in,
                                Vector<ValueType> &y_inout,
                                const ValueType alpha,
                                const ValueType beta) const;

private:
  /** \cond PRIVATE */
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
  typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
  typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

#ifdef WITH_TRILINOS
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType>> m_domainSpace;
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType>> m_rangeSpace;
#else
  unsigned int m_rowCount;
  unsigned int m_columnCount;
#endif

  blcluster *m_blockCluster;
  AhmedMblock **m_blocksL;
  AhmedMblock **m_blocksU;

  IndexPermutation m_domainPermutation;
  IndexPermutation m_rangePermutation;
  /** \endcond */
};

} // namespace Bempp

#endif
