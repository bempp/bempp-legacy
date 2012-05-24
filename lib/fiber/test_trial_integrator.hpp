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

#ifndef fiber_test_trial_integrator_hpp
#define fiber_test_trial_integrator_hpp

#include "scalar_traits.hpp"

#include <armadillo>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Basis;

/** \brief Integration of products of test and trial functions over elements. */
template <typename BasisFunctionType, typename ResultType>
class TestTrialIntegrator
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    virtual ~TestTrialIntegrator() {}

    virtual void integrate(
            const std::vector<int>& elementIndices,
            const Basis<BasisFunctionType>& testBasis,
            const Basis<BasisFunctionType>& trialBasis,
            arma::Cube<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
