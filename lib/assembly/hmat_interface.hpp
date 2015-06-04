// Copyright (C) 2011-2015 by the BEM++ Authors
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

#ifndef bempp_hmat_interface_hpp
#define bempp_hmat_interface_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../space/space.hpp"
#include "../hmat/geometry_interface.hpp"
#include "../hmat/block_cluster_tree.hpp"
#include "../fiber/scalar_traits.hpp"

#include <memory>
#include <vector>

namespace hmat {

/** \cond FORWARD_DECL */
class GeometryDataType;
}

namespace Bempp {

/** \brief This class provides between HMat and the grid.
 */
template <typename BasisFunctionType>
class SpaceHMatGeometryInterface : public hmat::GeometryInterface {
public:
  typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType
      CoordinateType;

  /** \brief Constructor */
  SpaceHMatGeometryInterface(const Space<BasisFunctionType> &space);

  /** \brief Obtain next element from the Geometry */
  shared_ptr<const hmat::GeometryDataType> next() override;

  /** \brief Number of geometric entities */
  std::size_t numberOfEntities() const override;
  void reset() override;

private:
  std::size_t m_counter;
  std::vector<BoundingBox<CoordinateType>> m_bemppBoundingBoxes;
};

/** \brief Generate a block cluster tree from a given pair of spaces. */
template <typename BasisFunctionType>
shared_ptr<hmat::DefaultBlockClusterTreeType>
generateBlockClusterTree(const Space<BasisFunctionType> &testSpace,
                         const Space<BasisFunctionType> &trialSpace,
                         const ParameterList &parameterList);
}

#endif
