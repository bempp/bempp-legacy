// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_numerical_test_trial_integrator_hpp
#define fiber_numerical_test_trial_integrator_hpp

#include "../common/common.hpp"

#include "test_trial_integrator.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType, typename ResultType> class TestTrialIntegral;
/** \endcond */

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
class NumericalTestTrialIntegrator :
        public TestTrialIntegrator<BasisFunctionType, ResultType>
{
public:
    typedef typename
    TestTrialIntegrator<BasisFunctionType, ResultType>::CoordinateType CoordinateType;

    NumericalTestTrialIntegrator(
            const arma::Mat<CoordinateType>& localQuadPoints,
            const std::vector<CoordinateType> quadWeights,
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const CollectionOfShapesetTransformations<CoordinateType>& testTransformations,
            const CollectionOfShapesetTransformations<CoordinateType>& trialTransformations,
            const TestTrialIntegral<BasisFunctionType, ResultType>& integral,
            const OpenClHandler& openClHandler);

    virtual void integrate(
            const std::vector<int>& elementIndices,
            const Shapeset<BasisFunctionType>& testShapeset,
            const Shapeset<BasisFunctionType>& trialShapeset,
            arma::Cube<ResultType>& result) const;

private:
    arma::Mat<CoordinateType> m_localQuadPoints;
    std::vector<CoordinateType> m_quadWeights;

    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;
    const CollectionOfShapesetTransformations<CoordinateType>& m_testTransformations;
    const CollectionOfShapesetTransformations<CoordinateType>& m_trialTransformations;
    const TestTrialIntegral<BasisFunctionType, ResultType>& m_integral;

    const OpenClHandler& m_openClHandler;
};

} // namespace Fiber

#include "numerical_test_trial_integrator_imp.hpp"

#endif
