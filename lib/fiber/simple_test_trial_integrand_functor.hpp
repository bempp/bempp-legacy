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

#ifndef fiber_simple_test_trial_integrand_functor_hpp
#define fiber_simple_test_trial_integrand_functor_hpp

#include "../common/common.hpp"

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"
#include "scalar_traits.hpp"

namespace Fiber
{

/** \ingroup functors
 *  \brief Functor evaluating the integrand of the standard identity operator.
 *
 *  This functor evaluates the integrand having the form of the scalar product
 *  of a single test basis function transformation and a single trial basis
 *  function transformation.
 */
template <typename BasisFunctionType_, typename ResultType_>
class SimpleTestTrialIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(size_t& geomDeps) const {
        // do nothing
    }

    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& /* geomData */,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& testValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& trialValues) const {
        // Assert that there is at least one test and trial transformation
        // and that the dimensions of the first pair agree
        assert(testValues.size() >= 1);
        assert(trialValues.size() >= 1);
        assert(testValues[0].extent(0) == trialValues[0].extent(0));

        const int transformationDim = testValues[0].extent(0);

        ResultType result = 0.;
        for (int dim = 0; dim < transformationDim; ++dim)
            result += conjugate(testValues[0](dim)) *
                    trialValues[0](dim);
        return result;
    }
};

} // namespace Fiber

#endif
