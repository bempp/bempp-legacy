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

#ifndef fiber_kernel_trial_integrator_hpp
#define fiber_kernel_trial_integrator_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"
#include "types.hpp"

#include <utility>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ValueType> class Shapeset;
/** \endcond */

/** \brief Integration over pairs of elements. */
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegrator {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
  typedef std::pair<int, int> PointElementIndexPair;

  virtual ~KernelTrialIntegrator() {}

  virtual void
  integrate(const std::vector<int> &pointIndices, int trialElementIndex,
            const Shapeset<BasisFunctionType> &trialShapeset,
            LocalDofIndex localTrialDofIndex,
            const std::vector<Matrix<ResultType> *> &result) const = 0;

  virtual void
  integrate(int pointIndex, int componentIndex,
            const std::vector<int> &trialElementIndices,
            const Shapeset<BasisFunctionType> &trialShapeset,
            const std::vector<Matrix<ResultType> *> &result) const = 0;

  virtual void
  integrate(const std::vector<PointElementIndexPair> &pointElementIndexPairs,
            const Shapeset<BasisFunctionType> &trialShapeset,
            const std::vector<Matrix<ResultType> *> &result) const = 0;
};

} // namespace Fiber

#endif
