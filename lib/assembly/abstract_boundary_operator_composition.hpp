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

#ifndef bempp_boundary_operator_composition_hpp
#define bempp_boundary_operator_composition_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator.hpp"
#include "boundary_operator.hpp"

#include "../common/shared_ptr.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \ingroup composite_boundary_operators
 *  \brief Composition of two abstract boundary operators.
 *
 *  This class represents a composition (product) of two boundary operators. */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperatorComposition :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    typedef typename Fiber::LocalAssemblerForOperators<ResultType>
    LocalAssembler;

    /** \brief Constructor.
     *
     *  Construct an operator representing the product \f$M \equiv L_1 L_2 : X
     *  \to Z\f$ of two boundary operators \f$L_1 : Y \to Z\f$ and \f$L_2 : X
     *  \to Y\f$.
     *
     *  \param[in] outer Operator \f$L_1\f$.
     *  \param[in] inner Operator \f$L_2\f$.
     *  \param[in] symmetry
     *    (Optional) Symmetry of the weak form of the composite operator.
     *    Can be any combination of the flags defined in the enumeration type
     *    Symmetry.
     *
     *  \note Both operators must be initialized and the range space of the
     *  operator \p inner must be identical with the domain space of the
     *  operator \p outer, otherwise an exception is thrown.
     *
     *  \todo Add a parameter setting the symmetry of the composite operator.
     */
    AbstractBoundaryOperatorComposition(
            const BoundaryOperator<BasisFunctionType, ResultType>& outer,
            const BoundaryOperator<BasisFunctionType, ResultType>& inner,
            int symmetry = NO_SYMMETRY);

    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleConversionOperator(const QuadratureStrategy& quadStrategy,
                               const AssemblyOptions& options);

private:
    BoundaryOperator<BasisFunctionType, ResultType> m_outer, m_inner;
};

} //namespace Bempp

#endif
