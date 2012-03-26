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

#ifndef fiber_separable_numerical_test_kernel_trial_integrator_hpp
#define fiber_separable_numerical_test_kernel_trial_integrator_hpp

#include "test_kernel_trial_integrator.hpp"
#include "raw_grid_geometry.hpp"
#include "opencl_handler.hpp"

namespace Fiber
{

template <typename ValueType, typename IndexType> class OpenClHandler;
template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename ValueType, typename GeometryFactory>
class SeparableNumericalTestKernelTrialIntegrator :
        public TestKernelTrialIntegrator<ValueType>
{
public:
    typedef typename TestKernelTrialIntegrator<ValueType>::ElementIndexPair
    ElementIndexPair;

    SeparableNumericalTestKernelTrialIntegrator(
            const arma::Mat<ValueType>& localTestQuadPoints,
            const arma::Mat<ValueType>& localTrialQuadPoints,
            const std::vector<ValueType> testQuadWeights,
            const std::vector<ValueType> trialQuadWeights,
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler);

    virtual void integrate(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            const Basis<ValueType>& basisA,
            const Basis<ValueType>& basisB,
            LocalDofIndex localDofIndexB,
            arma::Cube<ValueType>& result) const;

    virtual void integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis,
            arma::Cube<ValueType>& result) const;

private:
    virtual void integrateCpu(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            const Basis<ValueType>& basisA,
            const Basis<ValueType>& basisB,
            LocalDofIndex localDofIndexB,
            arma::Cube<ValueType>& result) const;

    virtual void integrateCl(
	    CallVariant callVariant,
	    const std::vector<int>& elementIndicesA,
	    int elementIndexB,
	    const Basis<ValueType>& basisA,
	    const Basis<ValueType>& basisB,
	    LocalDofIndex localDofIndexB,
	    arma::Cube<ValueType>& result) const;

    /**
     * \brief Returns an OpenCL code snippet containing the clIntegrate
     *   kernel function for integrating a single row or column
     */
    const std::string clStrIntegrateRowOrCol () const;

    arma::Mat<ValueType> m_localTestQuadPoints;
    arma::Mat<ValueType> m_localTrialQuadPoints;
    std::vector<ValueType> m_testQuadWeights;
    std::vector<ValueType> m_trialQuadWeights;

    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;

    const Expression<ValueType>& m_testExpression;
    const Kernel<ValueType>& m_kernel;
    const Expression<ValueType>& m_trialExpression;
    const OpenClHandler<ValueType,int>& m_openClHandler;
#ifdef WITH_OPENCL
    cl::Buffer *clTestQuadPoints;
    cl::Buffer *clTrialQuadPoints;
#endif
};

} // namespace Fiber

#include "separable_numerical_test_kernel_trial_integrator_imp.hpp"

#endif
