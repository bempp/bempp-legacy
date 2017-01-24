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

#ifndef fiber_hdiv_to_hcurl_function_value_functor_hpp
#define fiber_hdiv_to_hcurl_hcurl_function_value_functor_hpp

#include "../common/common.hpp"

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "shape_transformation_functor_wrappers.hpp"

namespace Fiber {

template <typename CoordinateType_>
class HdvivToHcurlFunctionValueElementaryFunctor {
public:
  typedef CoordinateType_ CoordinateType;

  int argumentDimension() const { return 2; }
  int resultDimension() const { return 3; }

  void addDependencies(size_t &basisDeps, size_t &geomDeps) const {
    basisDeps |= VALUES;
    geomDeps |= INTEGRATION_ELEMENTS | JACOBIANS_TRANSPOSED | NORMALS;
  }

  template <typename ValueType>
  void evaluate(const ConstBasisDataSlice<ValueType> &basisData,
                const ConstGeometricalDataSlice<CoordinateType> &geomData,
                _1dSliceOf3dArray<ValueType> &result) const {
    assert(basisData.componentCount() == argumentDimension());
    assert(result.extent(0) == resultDimension());
    boost::array<ValueType, 3> piola; // Value of function after Piola transform
    for (int rdim = 0; rdim < resultDimension(); ++rdim)
      piola[rdim] =
          (geomData.jacobianTransposed(0, rdim) * basisData.values(0) +
           geomData.jacobianTransposed(1, rdim) * basisData.values(1)) /
          geomData.integrationElement();
    // Now take the cross product with the normal direction

    result(0) = piola[1] * geomData.normal(2) - piola[2] * geomData.normal(1);
    result(1) = piola[2] * geomData.normal(0) - piola[0] * geomData.normal(2);
    result(2) = piola[0] * geomData.normal(1) - piola[1] * geomData.normal(0);
  }
};
// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
/** \ingroup functors
 *  \brief Functor calculating the value of a basis function from H(div). */
template <typename CoordinateType_>
class HdivToHcurlFunctionValueFunctor
    : public ElementaryShapeTransformationFunctorWrapper<
          HdivToHcurlFunctionValueElementaryFunctor<CoordinateType_>> {
public:
  typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
