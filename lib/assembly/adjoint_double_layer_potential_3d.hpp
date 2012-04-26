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

#ifndef bempp_adjoint_double_layer_potential_3d_hpp
#define bempp_adjoint_double_layer_potential_3d_hpp

#include "elementary_weakly_singular_integral_operator.hpp"
#include "../fiber/adjoint_double_layer_potential_3d_kernel.hpp"
#include "../fiber/scalar_function_value.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType = BasisFunctionType>
class AdjointDoubleLayerPotential3D :
        public ElementaryWeaklySingularIntegralOperator<BasisFunctionType, ResultType>
{
public:
    AdjointDoubleLayerPotential3D(const Space<BasisFunctionType>& testSpace,
                                  const Space<BasisFunctionType>& trialSpace);

private:
    virtual const Fiber::Kernel<ResultType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::Expression<BasisFunctionType>& testExpression() const {
        return m_expression;
    }

    virtual const Fiber::Expression<BasisFunctionType>& trialExpression() const {
        return m_expression;
    }

private:
    Fiber::AdjointDoubleLayerPotential3DKernel<ResultType> m_kernel;
    Fiber::ScalarFunctionValue<BasisFunctionType> m_expression;
};

} // namespace Bempp

#endif
