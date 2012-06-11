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

#include "test_trial_integrator.hpp"

namespace Fiber
{

class OpenClHandler;
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename CoordinateType> class RawGridGeometry;

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
            const CollectionOfBasisTransformations<CoordinateType>& testTransformations,
            const CollectionOfBasisTransformations<CoordinateType>& trialTransformations,
            const OpenClHandler& openClHandler);

    virtual void integrate(
            const std::vector<int>& elementIndices,
            const Basis<BasisFunctionType>& testBasis,
            const Basis<BasisFunctionType>& trialBasis,
            arma::Cube<ResultType>& result) const;

private:
    arma::Mat<CoordinateType> m_localQuadPoints;
    std::vector<CoordinateType> m_quadWeights;

    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;
    const CollectionOfBasisTransformations<CoordinateType>& m_testTransformations;
    const CollectionOfBasisTransformations<CoordinateType>& m_trialTransformations;

    const OpenClHandler& m_openClHandler;
};

} // namespace Fiber

#include "numerical_test_trial_integrator_imp.hpp"

#endif
