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

#ifndef fiber_standard_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_standard_local_assembler_for_integral_operators_on_surfaces_hpp

#include "local_assembler_for_operators.hpp"

#include "accuracy_options.hpp"
#include "array_2d.hpp"
#include "element_pair_topology.hpp"
#include "numerical_quadrature.hpp"
#include "parallelisation_options.hpp"
#include "shared_ptr.hpp"
#include "test_kernel_trial_integrator.hpp"

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <cstring>
#include <map>
#include <set>
#include <vector>

namespace Fiber
{

class OpenClHandler;

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
class StandardLocalAssemblerForIntegralOperatorsOnSurfaces :
        public LocalAssemblerForOperators<ResultType>
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    StandardLocalAssemblerForIntegralOperatorsOnSurfaces(
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Expression<CoordinateType> >& testExpression,
            const shared_ptr<const Kernel<KernelType> >& kernel,
            const shared_ptr<const Expression<CoordinateType> >& trialExpression,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals,
            const AccuracyOptions& accuracyOptions);
    virtual ~StandardLocalAssemblerForIntegralOperatorsOnSurfaces();

public:
    virtual void evaluateLocalWeakForms(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            LocalDofIndex localDofIndexB,
            std::vector<arma::Mat<ResultType> >& result);

    virtual void evaluateLocalWeakForms(
            const std::vector<int>& testElementIndices,
            const std::vector<int>& trialElementIndices,
            Fiber::Array2d<arma::Mat<ResultType> >& result);

    virtual void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Mat<ResultType> >& result);

private:
    typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> Integrator;
    typedef typename Integrator::ElementIndexPair ElementIndexPair;
    typedef std::set<ElementIndexPair> ElementIndexPairSet;
    typedef std::map<ElementIndexPair, arma::Mat<ResultType> > Cache;

    void checkConsistencyOfGeometryAndBases(
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& bases) const;

    bool testAndTrialGridsAreIdentical() const;

    void cacheSingularLocalWeakForms();
    void findPairsOfAdjacentElements(ElementIndexPairSet& pairs) const;
    void cacheLocalWeakForms(const ElementIndexPairSet& elementIndexPairs);

    const Integrator& selectIntegrator(
            int testElementIndex, int trialElementIndex);

    enum ElementType {
        TEST, TRIAL
    };

    int regularOrder(int elementIndex, ElementType elementType) const;
    int singularOrder(int elementIndex, ElementType elementType) const;

    const Integrator& getIntegrator(const DoubleQuadratureDescriptor& index);

private:
    typedef tbb::concurrent_unordered_map<DoubleQuadratureDescriptor,
    Integrator*> IntegratorMap;

private:
    shared_ptr<const GeometryFactory> m_testGeometryFactory;
    shared_ptr<const GeometryFactory> m_trialGeometryFactory;
    shared_ptr<const RawGridGeometry<CoordinateType> > m_testRawGeometry;
    shared_ptr<const RawGridGeometry<CoordinateType> > m_trialRawGeometry;
    shared_ptr<const std::vector<const Basis<BasisFunctionType>*> > m_testBases;
    shared_ptr<const std::vector<const Basis<BasisFunctionType>*> > m_trialBases;
    shared_ptr<const Expression<CoordinateType> > m_testExpression;
    shared_ptr<const Kernel<KernelType> > m_kernel;
    shared_ptr<const Expression<CoordinateType> > m_trialExpression;
    shared_ptr<const OpenClHandler> m_openClHandler;
    ParallelisationOptions m_parallelisationOptions;
    AccuracyOptions m_accuracyOptions;

    IntegratorMap m_TestKernelTrialIntegrators;
    Cache m_cache;
};

} // namespace Fiber

#endif
