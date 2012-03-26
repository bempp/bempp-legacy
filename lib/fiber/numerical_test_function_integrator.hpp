// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_numerical_test_function_integrator_hpp
#define fiber_numerical_test_function_integrator_hpp

#include "test_function_integrator.hpp"

namespace Fiber
{

template <typename ValueType, typename IndexType> class OpenClHandler;
template <typename ValueType> class Expression;
template <typename ValueType> class Function;
template <typename ValueType> class RawGridGeometry;

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename ValueType, typename GeometryFactory>
class NumericalTestFunctionIntegrator : public TestFunctionIntegrator<ValueType>
{
public:
    NumericalTestFunctionIntegrator(
            const arma::Mat<ValueType>& localQuadPoints,
            const std::vector<ValueType> quadWeights,
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const Expression<ValueType>& testExpression,
            const Function<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler);

    virtual void integrate(
            const std::vector<int>& elementIndices,
            const Basis<ValueType>& testBasis,
            arma::Mat<ValueType>& result) const;

private:
    arma::Mat<ValueType> m_localQuadPoints;
    std::vector<ValueType> m_quadWeights;

    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;
    const Expression<ValueType>& m_testExpression;
    const Function<ValueType>& m_function;

    const OpenClHandler<ValueType,int>& m_openClHandler;
};

} // namespace Fiber

#include "numerical_test_function_integrator_imp.hpp"

#endif
