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

#ifndef fiber_standard_evaluator_for_integral_operators_hpp
#define fiber_standard_evaluator_for_integral_operators_hpp

#include "evaluator_for_integral_operators.hpp"

#include <armadillo>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;
template <typename ValueType> class RawGridGeometry;
template <typename CoordinateType, typename IndexType> class OpenClHandler;

template <typename ValueType, typename GeometryFactory>
class StandardEvaluatorForIntegralOperators :
        public EvaluatorForIntegralOperators<ValueType>
{
public:
    typedef typename EvaluatorForIntegralOperators<ValueType>::Region Region;

    StandardEvaluatorForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const std::vector<std::vector<ValueType> >& argumentLocalCoefficients,
            ValueType multiplier,
            const OpenClHandler<ValueType, int>& openClHandler,
            const QuadratureOptions& quadratureOptions);

    virtual void evaluate(Region region,
                          const arma::Mat<ValueType>& points,
                          arma::Mat<ValueType>& result) const;

private:
    void cacheTrialData();
    void calcTrialData(
            Region region,
            int kernelTrialGeomDeps,
            Fiber::GeometricalData<ValueType>& trialGeomData,
            arma::Mat<ValueType>& weightedTrialExprValues) const;

    int quadOrder(const Fiber::Basis<ValueType>& basis, Region region) const;
    int farFieldQuadOrder(const Fiber::Basis<ValueType>& basis) const;
    int nearFieldQuadOrder(const Fiber::Basis<ValueType>& basis) const;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;
    const std::vector<const Basis<ValueType>*>& m_trialBases;
    const Kernel<ValueType>& m_kernel;
    const Expression<ValueType>& m_trialExpression;
    const std::vector<std::vector<ValueType> >& m_argumentLocalCoefficients;
    ValueType m_multiplier;
    const Fiber::OpenClHandler<ValueType,int>& m_openClHandler;
    const QuadratureOptions& m_quadratureOptions;

    Fiber::GeometricalData<ValueType> m_nearFieldTrialGeomData;
    Fiber::GeometricalData<ValueType> m_farFieldTrialGeomData;
    arma::Mat<ValueType> m_nearFieldWeightedTrialExprValues;
    arma::Mat<ValueType> m_farFieldWeightedTrialExprValues;
};

} // namespace Bempp

#include "standard_evaluator_for_integral_operators_imp.hpp"

#endif
