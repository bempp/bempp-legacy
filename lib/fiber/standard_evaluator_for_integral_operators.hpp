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

struct QuadratureOptions;
template <typename ValueType> class Basis;
template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;
template <typename CoordinateType> class RawGridGeometry;
template <typename CoordinateType, typename IndexType> class OpenClHandler;

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
class StandardEvaluatorForIntegralOperators :
        public EvaluatorForIntegralOperators<ResultType>
{
public:
    typedef EvaluatorForIntegralOperators<ResultType> Base;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::Region Region;

    StandardEvaluatorForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<KernelType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            const QuadratureOptions& quadratureOptions);

    virtual void evaluate(Region region,
                          const arma::Mat<CoordinateType>& points,
                          arma::Mat<ResultType>& result) const;

private:
    void cacheTrialData();
    void calcTrialData(
            Region region,
            int kernelTrialGeomDeps,
            Fiber::GeometricalData<CoordinateType>& trialGeomData,
            arma::Mat<ResultType>& weightedTrialExprValues) const;

    int quadOrder(const Fiber::Basis<BasisFunctionType>& basis, Region region) const;
    int farFieldQuadOrder(const Fiber::Basis<BasisFunctionType>& basis) const;
    int nearFieldQuadOrder(const Fiber::Basis<BasisFunctionType>& basis) const;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;
    const std::vector<const Basis<BasisFunctionType>*>& m_trialBases;
    const Kernel<KernelType>& m_kernel;
    const Expression<CoordinateType>& m_trialExpression;
    const std::vector<std::vector<ResultType> >& m_argumentLocalCoefficients;
    const Fiber::OpenClHandler<CoordinateType,int>& m_openClHandler;
    const QuadratureOptions& m_quadratureOptions;

    Fiber::GeometricalData<CoordinateType> m_nearFieldTrialGeomData;
    Fiber::GeometricalData<CoordinateType> m_farFieldTrialGeomData;
    arma::Mat<ResultType> m_nearFieldWeightedTrialExprValues;
    arma::Mat<ResultType> m_farFieldWeightedTrialExprValues;
};

} // namespace Bempp

#include "standard_evaluator_for_integral_operators_imp.hpp"

#endif
