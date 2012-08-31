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

#ifndef bempp_discrete_blocked_boundary_operator_hpp
#define bempp_discrete_blocked_boundary_operator_hpp

#include "bempp/common/config_trilinos.hpp"

#include "discrete_boundary_operator.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/_2d_array.hpp"

#ifdef WITH_TRILINOS
#include <Teuchos_RCP.hpp>
#endif // WITH_TRILINOS

namespace Bempp
{

/** \ingroup discrete_boundary_operators
 *  \brief Discrete boundary operator composed of multiple blocks stored separately.
 *
 *  This class represents a linear operator whose matrix
 *  \f[ L =
 *      \begin{bmatrix}
 *        L_{11} & L_{12} & \dots  & L_{1n} \\
 *        L_{21} & L_{22} & \dots  & L_{2n} \\
 *        \vdots & \vdots & \ddots & \vdots \\
 *        L_{m1} & L_{m2} & \dots  & L_{mn}
 *      \end{bmatrix}
 *  \f]
 *  is composed of \f$m \times n\f$ discrete boundary operators \f$L_{ij}\f$. */
template <typename ValueType>
class DiscreteBlockedBoundaryOperator : public DiscreteBoundaryOperator<ValueType>
{
public:
    typedef DiscreteBoundaryOperator<ValueType> Base;

    /** \brief Constructor.
     *
     *  \param[in] blocks
     *    2D array of shared pointers to the discrete boundary operators
     *    \f$L_{ij}\f$ making up the newly constructed blocked operator. A null
     *    shared pointer is equivalent to a discrete boundary operator with
     *    zero matrix.
     *  \param[in] rowCounts
     *    Vector whose <em>i</em>th element is the number of rows of the matrix
     *    of each operator in <em>i</em>th row of the array \p blocks.
     *  \param[in] columnCounts
     *    Vector whose <em>i</em>th element is the number of columns of the
     *    matrix of each operator in <em>i</em>th column of the array \p blocks.
     */
    DiscreteBlockedBoundaryOperator(
            const Fiber::_2dArray<shared_ptr<const Base> >& blocks,
            const std::vector<size_t>& rowCounts,
            const std::vector<size_t>& columnCounts);

    virtual unsigned int rowCount() const;
    virtual unsigned int columnCount() const;

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const;

    inline shared_ptr<const DiscreteBoundaryOperator<ValueType> > asDiscreteAcaBoundaryOperator(
                                                              double eps=1E-4,
                                                              int maximumRank=50) const{
        throw std::runtime_error("DiscreteBlockedBoundaryOperator::asDiscreteAcaBoundaryOperator:"
                                 "not implemented.");
    }


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
    Fiber::_2dArray<shared_ptr<const Base> > m_blocks;
    std::vector<size_t> m_rowCounts;
    std::vector<size_t> m_columnCounts;
#ifdef WITH_TRILINOS
    Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > m_rangeSpace;
#endif
    /** \endcond */
};

} // namespace Bempp

#endif
