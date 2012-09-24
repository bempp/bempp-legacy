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

#ifndef fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "local_assembler_for_operators.hpp"

#include "_2d_array.hpp"
#include "accuracy_options.hpp"
#include "element_pair_topology.hpp"
#include "numerical_quadrature.hpp"
#include "parallelization_options.hpp"
#include "shared_ptr.hpp"
#include "test_kernel_trial_integrator.hpp"
#include "verbosity_level.hpp"

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <cstring>
#include <map>
#include <set>
#include <vector>

namespace Fiber
{

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;
template <typename CoordinateType> class RawGridGeometry;
/** \endcond */

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
class DefaultLocalAssemblerForIntegralOperatorsOnSurfaces :
        public LocalAssemblerForOperators<ResultType>
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    DefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const CollectionOfKernels<KernelType> >& kernel,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType> >& integral,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            VerbosityLevel::Level verbosityLevel,
            bool cacheSingularIntegrals,
            const AccuracyOptionsEx& accuracyOptions);
    virtual ~DefaultLocalAssemblerForIntegralOperatorsOnSurfaces();

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
            Fiber::_2dArray<arma::Mat<ResultType> >& result);

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

    void getRegularOrders(int testElementIndex, int trialElementIndex,
                     int& testQuadOrder, int& trialQuadOrder) const;

    CoordinateType elementSizeSquared(
            int elementIndex, const RawGridGeometry<CoordinateType>& rawGeometry) const;
    arma::Col<CoordinateType> elementCenter(
            int elementIndex, const RawGridGeometry<CoordinateType>& rawGeometry) const;
    CoordinateType elementDistanceSquared(
            int testElementIndex, int trialElementIndex) const;

    void precalculateElementSizesAndCenters();
    void precalculateElementSizesAndCentersForSingleGrid(
            const RawGridGeometry<CoordinateType>& rawGeometry,
            std::vector<CoordinateType>& elementSizesSquared,
            arma::Mat<CoordinateType>& elementCenters);

private:
    typedef tbb::concurrent_unordered_map<DoubleQuadratureDescriptor,
    Integrator*> IntegratorMap;

private:
    /** \cond */
    shared_ptr<const GeometryFactory> m_testGeometryFactory;
    shared_ptr<const GeometryFactory> m_trialGeometryFactory;
    shared_ptr<const RawGridGeometry<CoordinateType> > m_testRawGeometry;
    shared_ptr<const RawGridGeometry<CoordinateType> > m_trialRawGeometry;
    shared_ptr<const std::vector<const Basis<BasisFunctionType>*> > m_testBases;
    shared_ptr<const std::vector<const Basis<BasisFunctionType>*> > m_trialBases;
    shared_ptr<const CollectionOfBasisTransformations<CoordinateType> > m_testTransformations;
    shared_ptr<const CollectionOfKernels<KernelType> > m_kernels;
    shared_ptr<const CollectionOfBasisTransformations<CoordinateType> > m_trialTransformations;
    shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType> > m_integral;
    shared_ptr<const OpenClHandler> m_openClHandler;
    ParallelizationOptions m_parallelizationOptions;
    VerbosityLevel::Level m_verbosityLevel;
    AccuracyOptionsEx m_accuracyOptions;

    IntegratorMap m_TestKernelTrialIntegrators;
    Cache m_cache;
    std::vector<CoordinateType> m_testElementSizesSquared;
    std::vector<CoordinateType> m_trialElementSizesSquared;
    arma::Mat<CoordinateType> m_testElementCenters;
    arma::Mat<CoordinateType> m_trialElementCenters;
};

} // namespace Fiber

#include "default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
