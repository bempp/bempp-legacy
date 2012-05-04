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

#ifndef bempp_dot_3d_double_layer_potential_hpp
#define bempp_dot_3d_double_layer_potential_hpp

#include "elementary_weakly_singular_integral_operator.hpp"
#include "../fiber/dot_3d_double_layer_potential_kernel.hpp"
#include "../fiber/scalar_function_value.hpp"
#include "../common/scalar_traits.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType = BasisFunctionType>
class Dot3dDoubleLayerPotential :
        public ElementaryWeaklySingularIntegralOperator<
        BasisFunctionType,
        ResultType, // typename ScalarTraits<ResultType>::RealType,
        ResultType>
{
  //typedef typename ScalarTraits<ResultType>::RealType KernelType;
    typedef ResultType KernelType;
    typedef ElementaryWeaklySingularIntegralOperator<
    BasisFunctionType, KernelType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    Dot3dDoubleLayerPotential(const Space<BasisFunctionType>& testSpace,
			      const Space<BasisFunctionType>& trialSpace,
			      KernelType waveNumber);

private:
    virtual const Fiber::Kernel<KernelType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::Expression<CoordinateType>& testExpression() const {
        return m_expression;
    }

    virtual const Fiber::Expression<CoordinateType>& trialExpression() const {
        return m_expression;
    }

private:
    Fiber::Dot3dDoubleLayerPotentialKernel<KernelType> m_kernel;
    Fiber::ScalarFunctionValue<CoordinateType> m_expression;
};

} // namespace Bempp

#endif
