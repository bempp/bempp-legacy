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

#ifndef fiber_standard_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_standard_local_assembler_for_integral_operators_on_surfaces_hpp

#include "local_assembler_for_integral_operators.hpp"

#include "array_2d.hpp"
#include "element_pair_topology.hpp"
#include "opencl_options.hpp"
#include "nonseparable_numerical_test_kernel_trial_integrator.hpp"
#include "numerical_quadrature.hpp"
#include "separable_numerical_test_kernel_trial_integrator.hpp"

#include <boost/static_assert.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <cstring>
#include <map>
#include <set>
#include <vector>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
class StandardLocalAssemblerForIntegralOperatorsOnSurfaces :
        public LocalAssemblerForIntegralOperators<ValueType>
{    
public:
    StandardLocalAssemblerForIntegralOperatorsOnSurfaces(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClOptions& openClOptions);

public:
    virtual void evaluateLocalWeakForms(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            LocalDofIndex localDofIndexB,
            std::vector<arma::Mat<ValueType> >& result);

    virtual void evaluateLocalWeakForms(
            const std::vector<int>& testElementIndices,
            const std::vector<int>& trialElementIndices,
            Fiber::Array2D<arma::Mat<ValueType> >& result);

private:
    const TestKernelTrialIntegrator<ValueType>& selectIntegrator(
            int testElementIndex, int trialElementIndex);

    int regularOrderIncrement(int elementIndex) const;
    int singularOrderIncrement(int elementIndex) const;

    const TestKernelTrialIntegrator<ValueType>& getIntegrator(
            const DoubleQuadratureDescriptor& index);

private:
    typedef boost::ptr_map<DoubleQuadratureDescriptor,
    TestKernelTrialIntegrator<ValueType> > IntegratorMap;

    typedef typename TestKernelTrialIntegrator<ValueType>::ElementIndexPair ElementIndexPair;
    typedef std::map<ElementIndexPair, arma::Mat<ValueType> > Cache;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;
    const std::vector<const Basis<ValueType>*>& m_testBases;
    const std::vector<const Basis<ValueType>*>& m_trialBases;
    const Expression<ValueType>& m_testExpression;
    const Kernel<ValueType>& m_kernel;
    const Expression<ValueType>& m_trialExpression;
    OpenClOptions m_openClOptions;

    IntegratorMap m_TestKernelTrialIntegrators;
    Cache m_cache;
};

} // namespace Fiber

#include "standard_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
