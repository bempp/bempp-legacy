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

#ifndef bempp_laplace_3d_double_layer_potential_hpp
#define bempp_laplace_3d_double_layer_potential_hpp

#include "elementary_singular_integral_operator.hpp"
#include "../fiber/expression_list.hpp"
#include "../fiber/laplace_3d_double_layer_potential_kernel.hpp"
#include "../fiber/scalar_function_value.hpp"
#include "../common/scalar_traits.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType = BasisFunctionType>
class Laplace3dDoubleLayerPotential :
        public ElementarySingularIntegralOperator<
        BasisFunctionType,
        typename ScalarTraits<ResultType>::RealType,
        ResultType>
{
    typedef typename ScalarTraits<ResultType>::RealType KernelType;
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType, KernelType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    Laplace3dDoubleLayerPotential(const Space<BasisFunctionType>& testSpace,
                                  const Space<BasisFunctionType>& trialSpace);

private:
    virtual const Fiber::Kernel<KernelType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::ExpressionList<ResultType>& testExpressionList() const {
        return m_expressionList;
    }

    virtual const Fiber::ExpressionList<ResultType>& trialExpressionList() const {
        return m_expressionList;
    }

private:
    Fiber::Laplace3dDoubleLayerPotentialKernel<KernelType> m_kernel;
    Fiber::ScalarFunctionValue<CoordinateType> m_expression;
    Fiber::ExpressionList<ResultType> m_expressionList;
};

} // namespace Bempp

#endif
