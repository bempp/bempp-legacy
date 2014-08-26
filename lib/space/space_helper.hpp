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

#ifndef bempp_space_utils_hpp
#define bempp_space_utils_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/types.hpp"

#include <vector>

namespace Bempp {

class GridView;
struct LocalDof;
template <typename CoordinateType> struct BoundingBox;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType> class SpaceHelper {
public:
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  static void getGlobalDofInterpolationPoints_defaultImplementation(
      const Space<BasisFunctionType> &space, arma::Mat<CoordinateType> &points);

  static void getNormalsAtGlobalDofInterpolationPoints_defaultImplementation(
      const Space<BasisFunctionType> &space,
      arma::Mat<CoordinateType> &normals);

  static void getGlobalDofBoundingBoxes_defaultImplementation(
      const GridView &view,
      const std::vector<std::vector<LocalDof>> &global2localDofs,
      std::vector<BoundingBox<CoordinateType>> &bboxes);

  static void getGlobalDofNormals_defaultImplementation(
      const GridView &view,
      const std::vector<std::vector<LocalDof>> &global2localDofs,
      std::vector<Point3D<CoordinateType>> &normals);

  static void initializeLocal2FlatLocalDofMap(
      size_t flatLocalDofCount,
      const std::vector<std::vector<GlobalDofIndex>> &local2globalDofs,
      std::vector<LocalDof> &flatLocal2localDofs);
};

} // namespace Bempp

#endif
