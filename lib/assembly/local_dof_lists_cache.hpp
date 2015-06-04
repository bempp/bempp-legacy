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

#ifndef bempp_local_dof_lists_cache_hpp
#define bempp_local_dof_lists_cache_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <vector>
#include <iostream>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *
 *  \brief Data used by WeakFormAcaAssemblyHelper to convert between
 *  H-matrix indices, global and local degrees of freedom. */
template <typename BasisFunctionType> struct LocalDofLists {
  /** \brief Type used to index matrices.
   *
   *  Equivalent to either GlobalDofIndex (if m_indexWithGlobalDofs is true)
   *  or FlatLocalDofIndex (if m_indexWithGlobalDofs is false). */
  typedef int DofIndex;
  std::vector<DofIndex> originalIndices;
  std::vector<int> elementIndices;
  std::vector<std::vector<LocalDofIndex>> localDofIndices;
  std::vector<std::vector<BasisFunctionType>> localDofWeights;
  std::vector<std::vector<int>> arrayIndices;
};

/** \ingroup weak_form_assembly_internal
 *
 *  \brief Cache of LocalDofLists objects. */
template <typename BasisFunctionType> class LocalDofListsCache {
public:
  LocalDofListsCache(const Space<BasisFunctionType> &space,
                     const std::vector<std::size_t> &p2o,
                     bool indexWithGlobalDofs);
  ~LocalDofListsCache();

  /** \brief Return the LocalDofLists object describing the DOFs corresponding
   * to
   *  AHMED matrix indices [start, start + indexCount). */
  shared_ptr<const LocalDofLists<BasisFunctionType>> get(int start,
                                                         int indexCount);

private:
  void findLocalDofs(
      int start, int indexCount,
      std::vector<typename LocalDofLists<BasisFunctionType>::DofIndex> &
          originalIndices,
      std::vector<int> &elementIndices,
      std::vector<std::vector<LocalDofIndex>> &localDofIndices,
      std::vector<std::vector<BasisFunctionType>> &localDofWeights,
      std::vector<std::vector<int>> &arrayIndices) const;

  void findLocalDofs(
      int index,
      std::vector<typename LocalDofLists<BasisFunctionType>::DofIndex> &
          originalIndices,
      std::vector<int> &elementIndices,
      std::vector<std::vector<LocalDofIndex>> &localDofIndices,
      std::vector<std::vector<BasisFunctionType>> &localDofWeights,
      std::vector<std::vector<int>> &arrayIndices) const;

private:
  /** \cond PRIVATE */
  const Space<BasisFunctionType> &m_space;
  const std::vector<std::size_t> &m_p2o;
  bool m_indexWithGlobalDofs;

  typedef tbb::concurrent_unordered_map<
      std::pair<int, int>, const LocalDofLists<BasisFunctionType> *>
      LocalDofListsMap;
  LocalDofListsMap m_map;
  /** \endcond */
};

} // namespace Bempp

#endif
