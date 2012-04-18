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

#ifndef fiber_test_kernel_trial_integrator_hpp
#define fiber_test_kernel_trial_integrator_hpp

#include "scalar_traits.hpp"
#include "types.hpp"

#include <armadillo>
#include <utility>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Basis;

/** \brief Integration over pairs of elements. */
template <typename BasisValueType, typename KernelValueType>
class TestKernelTrialIntegrator
{
public:
    typedef typename Coercion<BasisValueType, KernelValueType>::Type ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    typedef std::pair<int, int> ElementIndexPair;

    virtual ~TestKernelTrialIntegrator() {}

    virtual void integrate(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            const Basis<BasisValueType>& basisA,
            const Basis<BasisValueType>& basisB,
            LocalDofIndex localDofIndexB,
            arma::Cube<ResultType>& result) const = 0;

    virtual void integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<BasisValueType>& testBasis,
            const Basis<BasisValueType>& trialBasis,
            arma::Cube<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
