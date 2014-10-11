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

#include "laplace_3d_calderon_projector.hpp"

#include "../common/shared_ptr.hpp"

#include "laplace_3d_single_layer_boundary_operator.hpp"
#include "laplace_3d_double_layer_boundary_operator.hpp"
#include "laplace_3d_hypersingular_boundary_operator.hpp"
#include "identity_operator.hpp"
#include "blocked_operator_structure.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../assembly/context.hpp"

#include <boost/make_shared.hpp>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dExteriorCalderonProjector(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label) {

  typedef BoundaryOperator<BasisFunctionType, ResultType> BdOp;

  shared_ptr<const Space<BasisFunctionType>> internalHplusSpace =
      hplusSpace->discontinuousSpace(hplusSpace);

  BdOp internalSlp = laplace3dSingleLayerBoundaryOperator(
      context, internalHplusSpace, internalHplusSpace, internalHplusSpace,
      label + "_slp", SYMMETRIC);

  BdOp dlp = laplace3dDoubleLayerBoundaryOperator(context, hplusSpace,
                                                  hplusSpace, hminusSpace,
                                                  label + "_dlp", NO_SYMMETRY);

  BdOp adjDlp = adjoint(dlp);

  BdOp hyp = laplace3dHypersingularBoundaryOperator(
      context, hplusSpace, hminusSpace, hplusSpace, label + "_hyp", SYMMETRIC,
      internalSlp);

  BdOp idSpaceTransformation1 = identityOperator(
      context, hminusSpace, internalHplusSpace, internalHplusSpace);
  BdOp idSpaceTransformation2 =
      identityOperator(context, internalHplusSpace, hplusSpace, hminusSpace);
  BdOp idDouble =
      identityOperator(context, hplusSpace, hplusSpace, hminusSpace);
  BdOp idAdjDouble =
      identityOperator(context, hminusSpace, hminusSpace, hplusSpace);

  // Now Assemble the entries of the Calderon Projector

  BlockedOperatorStructure<BasisFunctionType, ResultType> structure;

  structure.setBlock(0, 0, .5 * idDouble + dlp);
  structure.setBlock(0, 1, -1. * idSpaceTransformation2 * internalSlp *
                               idSpaceTransformation1);
  structure.setBlock(1, 0, -1. * hyp);
  structure.setBlock(1, 1, .5 * idAdjDouble - adjDlp);

  return BlockedBoundaryOperator<BasisFunctionType, ResultType>(structure);
}

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dInteriorCalderonProjector(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label) {

  typedef BoundaryOperator<BasisFunctionType, ResultType> BdOp;

  shared_ptr<const Space<BasisFunctionType>> internalHplusSpace =
      hplusSpace->discontinuousSpace(hplusSpace);

  BdOp internalSlp = laplace3dSingleLayerBoundaryOperator(
      context, internalHplusSpace, internalHplusSpace, internalHplusSpace,
      label + "_slp", SYMMETRIC);

  BdOp dlp = laplace3dDoubleLayerBoundaryOperator(context, hplusSpace,
                                                  hplusSpace, hminusSpace,
                                                  label + "_dlp", NO_SYMMETRY);

  BdOp adjDlp = adjoint(dlp);

  BdOp hyp = laplace3dHypersingularBoundaryOperator(
      context, hplusSpace, hminusSpace, hplusSpace, label + "_hyp", SYMMETRIC,
      internalSlp);

  BdOp idSpaceTransformation1 = identityOperator(
      context, hminusSpace, internalHplusSpace, internalHplusSpace);
  BdOp idSpaceTransformation2 =
      identityOperator(context, internalHplusSpace, hplusSpace, hminusSpace);
  BdOp idDouble =
      identityOperator(context, hplusSpace, hplusSpace, hminusSpace);
  BdOp idAdjDouble =
      identityOperator(context, hminusSpace, hminusSpace, hplusSpace);

  // Now Assemble the entries of the Calderon Projector

  BlockedOperatorStructure<BasisFunctionType, ResultType> structure;

  structure.setBlock(0, 0, .5 * idDouble - dlp);
  structure.setBlock(0, 1, idSpaceTransformation2 * internalSlp *
                               idSpaceTransformation1);
  structure.setBlock(1, 0, hyp);
  structure.setBlock(1, 1, .5 * idAdjDouble + adjDlp);

  return BlockedBoundaryOperator<BasisFunctionType, ResultType>(structure);
}

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dExteriorCalderonProjector(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label) {

  shared_ptr<const Context<BasisFunctionType, ResultType>> context(
      new Context<BasisFunctionType, ResultType>(parameterList));
  return laplace3dExteriorCalderonProjector<BasisFunctionType, ResultType>(
      context, hminusSpace, hplusSpace, label);
}

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dInteriorCalderonProjector(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label) {

  shared_ptr<const Context<BasisFunctionType, ResultType>> context(
      new Context<BasisFunctionType, ResultType>(parameterList));
  return laplace3dInteriorCalderonProjector<BasisFunctionType, ResultType>(
      context, hminusSpace, hplusSpace, label);
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT)                       \
  template BlockedBoundaryOperator<BASIS, RESULT>                              \
  laplace3dExteriorCalderonProjector(                                          \
      const ParameterList &, const shared_ptr<const Space<BASIS>> &,           \
      const shared_ptr<const Space<BASIS>> &, const std::string &);            \
  template BlockedBoundaryOperator<BASIS, RESULT>                              \
  laplace3dInteriorCalderonProjector(                                          \
      const ParameterList &, const shared_ptr<const Space<BASIS>> &,           \
      const shared_ptr<const Space<BASIS>> &, const std::string &);            \
  template BlockedBoundaryOperator<BASIS, RESULT>                              \
  laplace3dExteriorCalderonProjector(                                          \
      const shared_ptr<const Context<BASIS, RESULT>> &,                        \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &, const std::string &);            \
  template BlockedBoundaryOperator<BASIS, RESULT>                              \
  laplace3dInteriorCalderonProjector(                                          \
      const shared_ptr<const Context<BASIS, RESULT>> &,                        \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &, const std::string &)

FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // Bempp
