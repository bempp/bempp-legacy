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

#include "bempp/common/config_trilinos.hpp"

#ifndef bempp_discrete_null_boundary_operator_hpp
#define bempp_discrete_null_boundary_operator_hpp

#include "../common/common.hpp"

#include "discrete_boundary_operator.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
#endif

namespace Bempp
{

/** \ingroup discrete_boundary_operators
 *  \brief Null (zero-valued) discrete boundary operator. */
template <typename ValueType>
class DiscreteNullBoundaryOperator :
        public DiscreteBoundaryOperator<ValueType>
{
public:
    /** \brief Constructor.
     *
     *  Construct a null discrete boundary operator of dimensions
     *  (\p rowCount_, \p columnCount_). */
    DiscreteNullBoundaryOperator(size_t rowCount_, size_t columnCount_);

    virtual void dump() const;

    virtual arma::Mat<ValueType> asMatrix() const;

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

#ifdef WITH_TRILINOS
public:
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;
#endif

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const;

private:
    /** \cond PRIVATE */
#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<ValueType> > m_rangeSpace;
#else
    size_t m_rowCount, m_columnCount;
#endif
    /** \endcond */
};

} // namespace Bempp

#endif
