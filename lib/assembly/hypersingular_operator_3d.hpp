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

#ifndef bempp_hypersingular_operator_3d_hpp
#define bempp_hypersingular_operator_3d_hpp

#include "elementary_weakly_singular_integral_operator.hpp"
#include "../fiber/single_layer_potential_3d_kernel.hpp"
#include "../fiber/surface_curl_3d.hpp"

namespace Bempp
{

// Hypersingular is obviously *not* weakly singular... but its
// weak form apparently can be integrated with methods devised
// for weakly singular operators.
// FIXME: maybe rename "weaklysingular" to "singular"
template <typename ValueType>
class HypersingularOperator3D :
        public ElementaryWeaklySingularIntegralOperator<ValueType>
{

public:
    HypersingularOperator3D(const Space<ValueType>& testSpace,
                            const Space<ValueType>& trialSpace);

private:
    virtual const Fiber::Kernel<ValueType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::Expression<ValueType>& testExpression() const {
        return m_expression;
    }

    virtual const Fiber::Expression<ValueType>& trialExpression() const {
        return m_expression;
    }

private:
    Fiber::SingleLayerPotential3DKernel<ValueType> m_kernel;
    Fiber::SurfaceCurl3D<ValueType> m_expression;
};

} // namespace Bempp

#endif
