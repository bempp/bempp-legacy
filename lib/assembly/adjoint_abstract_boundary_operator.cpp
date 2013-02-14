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

#include "adjoint_abstract_boundary_operator.hpp"

#include "context.hpp"
#include "transposed_discrete_boundary_operator.hpp"

#include "../common/complex_aux.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/type_traits/is_complex.hpp>
#include <stdexcept>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>::
AdjointAbstractBoundaryOperator(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        int symmetry) :
    Base(boundaryOp.dualToRange(),
         boundaryOp.domain() == boundaryOp.range() ?
             boundaryOp.dualToRange() :
             // assume that domain == dualToRange, we'll verify it
             // in the body of the constructor
             boundaryOp.range(),
         boundaryOp.domain(),
         "adj(" + boundaryOp.label() + ")",
         symmetry & AUTO_SYMMETRY ?
             boundaryOp.abstractOperator()->symmetry() :
             symmetry),
    m_operator(boundaryOp)
{
    if (boost::is_complex<BasisFunctionType>())
        throw std::logic_error(
            "AdjointAbstractBoundaryOperator(): Taking the adjoint of operators "
            "acting on complex-valued basis functions is not supported");
    if (boundaryOp.domain() != boundaryOp.range() &&
            boundaryOp.domain() != boundaryOp.dualToRange())
        throw std::runtime_error(
                "AdjointAbstractBoundaryOperator::"
                "AdjointAbstractBoundaryOperator(): "
                "Dual to the domain of the operator to invert cannot be determined "
                "since the domain is different from both "
                "the range and the space dual to its range");
}

template <typename BasisFunctionType, typename ResultType>
AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>::
AdjointAbstractBoundaryOperator(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
	const shared_ptr<const Space<BasisFunctionType> >& range,
        int symmetry) :
    Base(boundaryOp.dualToRange(),
	 range,
         boundaryOp.domain(),
         "adj(" + boundaryOp.label() + ")",
         symmetry & AUTO_SYMMETRY ?
             boundaryOp.abstractOperator()->symmetry() :
             symmetry),
    m_operator(boundaryOp)
{
    if (boost::is_complex<BasisFunctionType>())
        throw std::logic_error(
            "AdjointAbstractBoundaryOperator(): Taking the adjoint of operators "
            "acting on complex-valued basis functions is not supported");
}


template <typename BasisFunctionType, typename ResultType>
bool AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>::isLocal() const
{
    return m_operator.abstractOperator()->isLocal();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const Context<BasisFunctionType, ResultType>& context) const
{
    return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                new TransposedDiscreteBoundaryOperator<ResultType>(
                    TRANSPOSE, m_operator.weakForm()));
}

// #define FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_REAL_BASIS_ONLY(CLASSNAME) \
//     FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_REAL(CLASSNAME); \
//     FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_SP_REAL_COMPLEX(CLASSNAME); \
//     FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_REAL(CLASSNAME); \
//     FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_DP_REAL_COMPLEX(CLASSNAME)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AdjointAbstractBoundaryOperator);

} // namespace Bempp
