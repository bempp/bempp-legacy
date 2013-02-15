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

#ifndef fiber_local_assembler_for_potential_operators_hpp
#define fiber_local_assembler_for_potential_operators_hpp

#include "../common/common.hpp"

#include "_3d_array.hpp"
#include "scalar_traits.hpp"
#include "types.hpp"

#include <vector>

namespace Fiber
{

const int ALL_COMPONENTS = -1;

/** \brief Abstract interface of a local assembler for potential operators.

  A local assembler constructs local (element-by-element) discrete
  representations of potential operators, i.e. small arrays such as

  \f[
  P_{ij} = \int_{E_j} K(x_i, y) f_j(y) \, \mathrm{d}E_j,
  \]

  where \f$\{x_i\}\f$ is a set of points, \f$E_j\f$ a set of elements and \f$f_j\f$
  a set of functions living on the corresponding elements.

  Concrete subclasses of this interface class automatically select integrators
  well-adapted for evaluation of specific integrals occurring in matrix
  representations of different classes of operators.
 */
template <typename ResultType>
class LocalAssemblerForPotentialOperators
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    virtual ~LocalAssemblerForPotentialOperators() {}

    /** \brief Assemble local contributions.

    If \p nominalDistance is nonnegative, it is taken as the distance between
    all point-element pairs for the purposes of selecting the quadrature method.
    Otherwise the interelement distance is calculated separately for each
    point-element pair. */
    virtual void evaluateLocalContributions(
            const std::vector<int>& pointIndices,
            int trialElementIndex,
            LocalDofIndex localTrialDofIndex,
            std::vector<arma::Mat<ResultType> >& result,
            CoordinateType nominalDistance = -1.) = 0;

    /** \brief Assemble local contributions.

    If \p nominalDistance is nonnegative, it is taken as the distance between
    all point-element pairs for the purposes of selecting the quadrature method.
    Otherwise the interelement distance is calculated separately for each
    point-element pair. */
    virtual void evaluateLocalContributions(
            int pointIndex,
            int componentIndex,
            const std::vector<int>& trialElementIndices,
            std::vector<arma::Mat<ResultType> >& result,
            CoordinateType nominalDistance = -1.) = 0;

    /** \brief Assemble local contributions.

    If \p nominalDistance is nonnegative, it is taken as the distance between
    all point-element pairs for the purposes of selecting the quadrature method.
    Otherwise the interelement distance is calculated separately for each
    point-element pair. */
    virtual void evaluateLocalContributions(
            const std::vector<int>& pointIndices,
            const std::vector<int>& trialElementIndices,
            Fiber::_2dArray<arma::Mat<ResultType> >& result,
            CoordinateType nominalDistance = -1.) = 0;

    virtual int resultDimension() const = 0;

    virtual CoordinateType estimateRelativeScale(CoordinateType minDist) const = 0;
};

} // namespace Fiber

#endif
