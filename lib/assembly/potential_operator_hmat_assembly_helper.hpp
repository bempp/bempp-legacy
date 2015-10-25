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

#ifndef bempp_potential_operator_hmat_assembly_helper_hpp
#define bempp_potential_operator_hmat_assembly_helper_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../common/eigen_support.hpp"
#include "../hmat/hmatrix.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../hmat/block_cluster_tree.hpp"

#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForPotentialOperators;
/** \endcond */
}

namespace Bempp {

/** \cond FORWARD_DECL */
class ComponentListsCache;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class LocalDofListsCache;
template <typename BasisFunctionType> class Space;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief Assembly helper called by HMat for Potential Operators.
 */

template <typename BasisFunctionType, typename ResultType>
class PotentialOperatorHMatAssemblyHelper
    : public hmat::DataAccessor<ResultType, 2> {
public:
  typedef Fiber::LocalAssemblerForPotentialOperators<ResultType> LocalAssembler;
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
  typedef CoordinateType MagnitudeType;

  PotentialOperatorHMatAssemblyHelper(
      const Matrix<CoordinateType> &points,
      const Space<BasisFunctionType> &trialSpace,
      const shared_ptr<const hmat::DefaultBlockClusterTreeType>
          &blockClusterTree,
      LocalAssembler &assembler,
      const ParameterList& parameterList);

  void computeMatrixBlock(
      const hmat::IndexRangeType &testIndexRange,
      const hmat::IndexRangeType &trialIndexRange,
      const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
      Matrix<ResultType> &data) const override;

  double
  scale(const hmat::DefaultBlockClusterTreeNodeType &node) const override;

private:
  MagnitudeType estimateMinimumDistance(
      const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode) const;

  const Matrix<CoordinateType> &m_points;
  const Space<BasisFunctionType> &m_trialSpace;
  const shared_ptr<const hmat::DefaultBlockClusterTreeType> m_blockClusterTree;
  LocalAssembler &m_assembler;
  const ParameterList& m_parameterList;
  int m_componentCount;

  shared_ptr<LocalDofListsCache<BasisFunctionType>> m_trialDofListsCache;
  shared_ptr<ComponentListsCache> m_componentListsCache;

};

} // namespace Bempp

#endif
