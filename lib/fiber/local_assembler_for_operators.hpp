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

#ifndef fiber_local_assembler_for_operators_hpp
#define fiber_local_assembler_for_operators_hpp

#include "../common/common.hpp"

#include "_2d_array.hpp"
#include "scalar_traits.hpp"
#include "types.hpp"

#include <vector>

namespace Fiber
{

/** \brief Abstract interface of a local assembler for operators.

  A local assembler constructs local (element-by-element) matrix
  representations of linear operators used in boundary-element methods.

  Typically, a separate local assembler must be constructed for each elementary
  integral operator and for the identity operator.

  Concrete subclasses of this interface class automatically select integrators
  well-adapted for evaluation of specific integrals occurring in matrix
  representations of different classes of operators.
 */
template <typename ResultType>
class LocalAssemblerForOperators
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    virtual ~LocalAssemblerForOperators() {}

//    /** \brief Number of dimensions over which integration takes place. */
//    int elementDimension() const {
//        asImp().elementDimension();
//    }

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
    virtual void evaluateLocalWeakForms(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            LocalDofIndex localDofIndexB,
            std::vector<arma::Mat<ResultType> >& result,
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
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& testElementIndices,
            const std::vector<int>& trialElementIndices,
            Fiber::_2dArray<arma::Mat<ResultType> >& result,
            CoordinateType nominalDistance = -1.) = 0;

    /** \brief Assemble local weak forms.

      \param[in]  elementIndices Vector of element indices.
      \param[out] result         Vector of weak forms of the operator on
                                 element pairs
                                 (\p element(\p i), \p element(\p i))
                                 for \p i in \p elementIndices. */
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Mat<ResultType> >& result) = 0;

    virtual CoordinateType estimateRelativeScale(CoordinateType minDist) const = 0;
};

} // namespace Fiber

#endif
