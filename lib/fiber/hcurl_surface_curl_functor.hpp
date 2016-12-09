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

#ifndef fiber_hcurl_surface_curl_functor_hpp
#define fiber_hcurl_surface_curl_functor_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "shape_transformation_functor_wrappers.hpp"

namespace Fiber {

template <typename CoordinateType_> class ElementaryHcurlSurfaceCurlFunctor {
public:
  typedef CoordinateType_ CoordinateType;

  int argumentDimension() const { return 2; }
  int resultDimension() const { return 1; }

  void addDependencies(size_t &basisDeps, size_t &geomDeps) const {
    basisDeps |= DERIVATIVES;
    geomDeps |= JACOBIAN_INVERSES_TRANSPOSED | NORMALS;
  }

  template <typename ValueType>
  void evaluate(const ConstBasisDataSlice<ValueType> &basisData,
                const ConstGeometricalDataSlice<CoordinateType> &geomData,
                _1dSliceOf3dArray<ValueType> &result) const {
    assert(basisData.componentCount() == argumentDimension());
    assert(result.extent(0) == resultDimension());

    const int cdim = 3;
    const int mydim = 2;

    Eigen::Matrix<CoordinateType, cdim, mydim> jacobianInverseTransposed;
    Eigen::Matrix<ValueType, mydim, mydim> referenceDerivative;
    Eigen::Matrix<ValueType, cdim, cdim> derivative;
    Eigen::Matrix<ValueType, cdim, 1> curlValue;

    for (int i = 0; i < cdim; ++i)
      for (int j = 0; j < mydim; ++j)
        jacobianInverseTransposed(i, j) =
            geomData.jacobianInverseTransposed(i, j);

    for (int i = 0; i < mydim; ++i)
      for (int j = 0; j < mydim; ++j)
        referenceDerivative(i, j) = basisData.derivatives(i, j);

    derivative = jacobianInverseTransposed * referenceDerivative *
                 jacobianInverseTransposed.transpose();

    curlValue(0) = derivative(2, 1) - derivative(1, 2);
    curlValue(1) = derivative(0, 2) - derivative(2, 0);
    curlValue(2) = derivative(1, 0) - derivative(0, 1);

    result(0) = curlValue(0) * geomData.normal(0) +
                curlValue(1) * geomData.normal(1) +
                curlValue(2) * geomData.normal(2);
  }
};
// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
/** \ingroup functors
 *  \brief Functor calculating the value of a basis function from H(div). */
template <typename CoordinateType_>
class HcurlSurfaceCurlFunctor
    : public ElementaryShapeTransformationFunctorWrapper<
          ElementaryHcurlSurfaceCurlFunctor<CoordinateType_>> {
public:
  typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
