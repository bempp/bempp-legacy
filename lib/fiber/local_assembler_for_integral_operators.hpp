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

#ifndef fiber_local_assembler_for_integral_operators_hpp
#define fiber_local_assembler_for_integral_operators_hpp

#include "../common/common.hpp"

#include "_2d_array.hpp"
#include "scalar_traits.hpp"
#include "types.hpp"

#include <vector>

namespace Fiber {

/** \brief Abstract interface of a local assembler for integral operators.

  Let \f$A\f$ be an integral boundary operator. A local assembler associated
  with this operator constructs matrices \f$A_{mn}\f$, where \f$m\f$ is a test
  element number and \f$n\f$ a trial element number, with entries
  \f$(A_{mn})_{ij} = (u_{mi}, A v_{nj}\f$) where \f$(\cdot, \cdot)\f$ is the
  standard scalar product, \f$(u_{mi})_{i=1}^m\f$ are the element-level test
  functions defined on element \f$m\f$ and \f$(v_{nj})_{j=1}^n\f$ are the
  element-level trial functions defined on element \f$n\f$.

  The local assembler is responsible for choosing an appropriate way of
  evaluating
  the necessary integrals. */
template <typename ResultType> class LocalAssemblerForIntegralOperators {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  virtual ~LocalAssemblerForIntegralOperators() {}

  /** \brief Assemble local weak forms.

  In this overload, a "column" of local weak forms is assembled. More
  specifically, on exit \p result is a vector of local weak forms corresponding
  to the following pairs of elements:

  - if \p callVariant is \p TEST_TRIAL, all pairs (\p elementA, \p elementB)
  for \p elementA in \p elementsA;

  - if \p callVariant is \p TRIAL_TEST, all pairs (\p elementB, \p elementA)
  for \p elementA in \p elementsA.

  Unless \p localDofIndexB is set to \p ALL_DOFS, only entries corresponding
  to the (\p localDofIndexB)th local DOF on \p elementB are calculated.

  If \p nominalDistance is nonnegative, it is taken as the distance between
  all element pairs for the purposes of selecting the quadrature method.
  Otherwise the interelement distance is calculated separately for each
  element pair. */
  virtual void
  evaluateLocalWeakForms(CallVariant callVariant,
                         const std::vector<int> &elementIndicesA,
                         int elementIndexB, LocalDofIndex localDofIndexB,
                         std::vector<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.) = 0;

  /** \brief Assemble local weak forms.

  This overload constructs and assigns to the output parameter \p result the
  2D array of local weak forms corresponding to all pairs (\p testElement, \p
  trialElement) with \p testElement in \p testElementIndices and \p
  trialElement in \p trialElementIndices.

  This function should be used primarily for small blocks of elements lying
  close to each other.

  If \p nominalDistance is nonnegative, it is taken as the distance between
  all element pairs for the purposes of selecting the quadrature method.
  Otherwise the interelement distance is calculated separately for each
  element pair. */
  virtual void
  evaluateLocalWeakForms(const std::vector<int> &testElementIndices,
                         const std::vector<int> &trialElementIndices,
                         Fiber::_2dArray<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.) = 0;

  /** \brief Estimate how fast the entries in the matrix of
   *  this operator decay with interelement distance.
   *
   *  This function should return a rough (!) estimate of the magnitude of
   *  the entries in the matrix of this operator associated with a pair of
   *  test and trial elements lying a distance minDist from each other,
   *  relative to the estimated magnitude of entries associated with
   *  coincident elements. If no good estimate can be given, the function can
   *  return 1.
   *
   *  \note This function is called during ACA to detect block clusters on
   *  which the operator becomes very small and so can be safely approximated
   *  with 0. */
  virtual CoordinateType
  estimateRelativeScale(CoordinateType minDist) const = 0;
};

} // namespace Fiber

#endif
