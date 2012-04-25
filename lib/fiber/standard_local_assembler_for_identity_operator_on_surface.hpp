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

#ifndef fiber_standard_local_assembler_for_identity_operator_on_surface_hpp
#define fiber_standard_local_assembler_for_identity_operator_on_surface_hpp

#include "local_assembler_for_operators.hpp"
#include "numerical_quadrature.hpp"
#include "numerical_test_trial_integrator.hpp"

#include <armadillo>
#include <boost/static_assert.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace Fiber
{

template <typename CoordinateType, typename IndexType> class OpenClHandler;

template <typename BasisValueType, typename ResultType, typename GeometryFactory>
class StandardLocalAssemblerForIdentityOperatorOnSurface :
    public LocalAssemblerForOperators<ResultType>
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    StandardLocalAssemblerForIdentityOperatorOnSurface(
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const std::vector<const Basis<BasisValueType>*>& testBases,
        const std::vector<const Basis<BasisValueType>*>& trialBases,
        const Expression<BasisValueType>& testExpression,
        const Expression<BasisValueType>& trialExpression,
        const OpenClHandler<CoordinateType, int>& openClHandler);

    virtual void evaluateLocalWeakForms(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        LocalDofIndex localDofIndexB,
        std::vector<arma::Mat<ResultType> >& result);

    virtual void evaluateLocalWeakForms(
        const std::vector<int>& testElementIndices,
        const std::vector<int>& trialElementIndices,
        Fiber::Array2D<arma::Mat<ResultType> >& result);

    virtual void evaluateLocalWeakForms(
        const std::vector<int>& elementIndices,
        std::vector<arma::Mat<ResultType> >& result);

private:
    const TestTrialIntegrator<BasisValueType, ResultType>&
    selectIntegrator(int elementIndex);

    const TestTrialIntegrator<BasisValueType, ResultType>& getIntegrator(
        const SingleQuadratureDescriptor& desc);
private:
    typedef boost::ptr_map<SingleQuadratureDescriptor,
            TestTrialIntegrator<BasisValueType, ResultType> > IntegratorMap;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;
    const std::vector<const Basis<BasisValueType>*>& m_testBases;
    const std::vector<const Basis<BasisValueType>*>& m_trialBases;
    const Expression<BasisValueType>& m_testExpression;
    const Expression<BasisValueType>& m_trialExpression;
    const OpenClHandler<CoordinateType, int>& m_openClHandler;

    IntegratorMap m_testTrialIntegrators;
};

} // namespace Fiber

#include "standard_local_assembler_for_identity_operator_on_surface_imp.hpp"

#endif
