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

#include "helmholtz_3d_single_layer_boundary_operator.hpp"

#include "modified_helmholtz_3d_single_layer_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
                 typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dSingleLayerBoundaryOperator(
    const shared_ptr<const Context<
        BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label, int symmetry, bool useInterpolation,
    int interpPtsPerWavelength) {
  typedef typename ScalarTraits<BasisFunctionType>::ComplexType ComplexType;
  return modifiedHelmholtz3dSingleLayerBoundaryOperator<
      BasisFunctionType, ComplexType, ComplexType>(
      context, domain, range, dualToRange, waveNumber / ComplexType(0., 1.),
      label, symmetry, useInterpolation, interpPtsPerWavelength);
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS)                               \
  template BoundaryOperator<BASIS, ScalarTraits<BASIS>::ComplexType>           \
  helmholtz3dSingleLayerBoundaryOperator(                                      \
      const shared_ptr<                                                        \
          const Context<BASIS, ScalarTraits<BASIS>::ComplexType>> &,           \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &,                                  \
      ScalarTraits<BASIS>::ComplexType, const std::string &, int, bool, int)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // namespace Bempp
