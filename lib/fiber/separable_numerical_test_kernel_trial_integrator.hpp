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

#include "config_opencl.hpp"

#ifndef fiber_separable_numerical_test_kernel_trial_integrator_hpp
#define fiber_separable_numerical_test_kernel_trial_integrator_hpp

#include "test_kernel_trial_integrator.hpp"

namespace Fiber
{

template <typename ValueType, typename IndexType> class OpenClHandler;
template <typename CoordinateType> class Expression;
template <typename ValueType> class Kernel;
template <typename CoordinateType> class RawGridGeometry;

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
class SeparableNumericalTestKernelTrialIntegrator :
        public TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
{
public:
    typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> Base;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::ElementIndexPair ElementIndexPair;

    SeparableNumericalTestKernelTrialIntegrator(
            const arma::Mat<CoordinateType>& localTestQuadPoints,
            const arma::Mat<CoordinateType>& localTrialQuadPoints,
            const std::vector<CoordinateType> testQuadWeights,
            const std::vector<CoordinateType> trialQuadWeights,
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const Expression<CoordinateType>& testExpression,
            const Kernel<KernelType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler);

    virtual ~SeparableNumericalTestKernelTrialIntegrator ();

    virtual void integrate(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            const Basis<BasisFunctionType>& basisA,
            const Basis<BasisFunctionType>& basisB,
            LocalDofIndex localDofIndexB,
            arma::Cube<ResultType>& result) const;

    virtual void integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<BasisFunctionType>& testBasis,
            const Basis<BasisFunctionType>& trialBasis,
            arma::Cube<ResultType>& result) const;

private:
    void integrateCpu(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            const Basis<BasisFunctionType>& basisA,
            const Basis<BasisFunctionType>& basisB,
            LocalDofIndex localDofIndexB,
            arma::Cube<ResultType>& result) const;

    void integrateCl(
	    CallVariant callVariant,
	    const std::vector<int>& elementIndicesA,
	    int elementIndexB,
        const Basis<BasisFunctionType>& basisA,
        const Basis<BasisFunctionType>& basisB,
	    LocalDofIndex localDofIndexB,
        arma::Cube<ResultType>& result) const;

    void integrateCpu(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<BasisFunctionType>& testBasis,
            const Basis<BasisFunctionType>& trialBasis,
            arma::Cube<ResultType>& result) const;

    void integrateCl(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<BasisFunctionType>& testBasis,
            const Basis<BasisFunctionType>& trialBasis,
            arma::Cube<ResultType>& result) const;

    /**
     * \brief Returns an OpenCL code snippet containing the clIntegrate
     *   kernel function for integrating a single row or column
     */
    const std::pair<const char*, int> clStrIntegrateRowOrCol () const;

    arma::Mat<CoordinateType> m_localTestQuadPoints;
    arma::Mat<CoordinateType> m_localTrialQuadPoints;
    std::vector<CoordinateType> m_testQuadWeights;
    std::vector<CoordinateType> m_trialQuadWeights;

    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;

    const Expression<CoordinateType>& m_testExpression;
    const Kernel<KernelType>& m_kernel;
    const Expression<CoordinateType>& m_trialExpression;
    const OpenClHandler<CoordinateType, int>& m_openClHandler;
#ifdef WITH_OPENCL
    cl::Buffer *clTestQuadPoints;
    cl::Buffer *clTrialQuadPoints;
    cl::Buffer *clTestQuadWeights;
    cl::Buffer *clTrialQuadWeights;
#endif
};

} // namespace Fiber

#include "separable_numerical_test_kernel_trial_integrator_imp.hpp"

#endif
