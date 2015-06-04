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

#ifndef bempp_cluster_construction_helper_hpp
#define bempp_cluster_construction_helper_hpp

#include "../common/common.hpp"

#include "../common/eigen_support.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"
#include "ahmed_aux_fwd.hpp"

#include <memory>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
template <typename CoordinateType> struct AhmedDofWrapper;
template <typename CoordinateType> struct Point3D;
template <typename CoordinateType> struct BoundingBox;
class AcaOptions;
class IndexPermutation;
/** \endcond */

template <typename BasisFunctionType> class ClusterConstructionHelper {
private:
  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
  typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
  typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;
  typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;

  static void constructBemCluster(const Space<BasisFunctionType> &space,
                                  bool indexWithGlobalDofs,
                                  const AcaOptions &acaOptions,
                                  shared_ptr<AhmedBemCluster> &cluster,
                                  shared_ptr<IndexPermutation> &o2p,
                                  shared_ptr<IndexPermutation> &p2o);

  static void constructBemCluster(const Matrix<CoordinateType> &points,
                                  int componentCount,
                                  const AcaOptions &acaOptions,
                                  shared_ptr<AhmedBemCluster> &cluster,
                                  shared_ptr<IndexPermutation> &o2p,
                                  shared_ptr<IndexPermutation> &p2o);

  static std::unique_ptr<AhmedBemBlcluster> constructBemBlockCluster(
      const AcaOptions &acaOptions, bool symmetric,
      AhmedBemCluster &testCluster, AhmedBemCluster &trialCluster,
      bool useStrongAdmissibilityCondition, unsigned int &blockCount);

  static void truncateBemBlockCluster(blcluster *cluster,
                                      const blcluster *refCluster);

  static void
  getComponentDofPositions(const Matrix<CoordinateType> &points,
                           int componentCount,
                           std::vector<Point3D<CoordinateType>> &positions);

  static void getComponentBoundingBoxes(
      const Matrix<CoordinateType> &points, int componentCount,
      std::vector<BoundingBox<CoordinateType>> &boundingBoxes);
};

} // namespace Bempp

#endif
