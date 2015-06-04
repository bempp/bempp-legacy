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

#include "scaled_abstract_boundary_operator.hpp"

#include "context.hpp"
#include "scaled_discrete_boundary_operator.hpp"

#include "../common/complex_aux.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

namespace {

template <typename ResultType>
int autoSymmetry(int origSymmetry, ResultType weight) {
  if (origSymmetry & HERMITIAN)
    if (imagPart(weight) != 0.)
      return origSymmetry & ~HERMITIAN;
  return origSymmetry;
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType>::
    ScaledAbstractBoundaryOperator(
        ResultType multiplier_,
        const BoundaryOperator<BasisFunctionType, ResultType> &multiplicand_,
        int symmetry)
    : Base(multiplicand_.domain(), multiplicand_.range(),
           multiplicand_.dualToRange(),
           "" + toString(multiplier_) + " * (" + multiplicand_.label() + ")",
           symmetry & AUTO_SYMMETRY
               ? autoSymmetry(throwIfUninitialized(
                                  multiplicand_,
                                  "ScaledAbstractBoundaryOperator::"
                                  "ScaledAbstractBoundaryOperator(): "
                                  "the boundary operator to be scaled must "
                                  "be initialized")
                                  .abstractOperator()
                                  ->symmetry(),
                              multiplier_)
               : symmetry),
      m_multiplier(multiplier_), m_multiplicand(multiplicand_) {}

template <typename BasisFunctionType, typename ResultType>
bool ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType>::isLocal()
    const {
  return m_multiplicand.abstractOperator()->isLocal();
}

template <typename BasisFunctionType_, typename ResultType_>
ResultType_ ScaledAbstractBoundaryOperator<BasisFunctionType_,
                                           ResultType_>::multiplier() const {
  return m_multiplier;
}

template <typename BasisFunctionType_, typename ResultType_>
BoundaryOperator<BasisFunctionType_, ResultType_>
ScaledAbstractBoundaryOperator<BasisFunctionType_, ResultType_>::multiplicand()
    const {
  return m_multiplicand;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    ScaledAbstractBoundaryOperator);

} // namespace Bempp
