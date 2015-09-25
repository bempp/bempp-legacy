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

#ifndef bempp_component_lists_cache_hpp
#define bempp_component_lists_cache_hpp

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

/** \ingroup potential_assembly_internal
 *
 *  \brief Data used by PotentialOperatorAcaAssemblyHelper to convert between
 *  H-matrix indices, point and component indices. */
struct ComponentLists {
  std::vector<int> pointIndices;
  std::vector<std::vector<int>> componentIndices;
  std::vector<std::vector<int>> arrayIndices;
};

/** \ingroup potential_assembly_internal
 *
 *  \brief Cache of ComponentLists objects. */
class ComponentListsCache {
public:
  ComponentListsCache(const std::vector<std::size_t> &p2o, int componentCount);
  ~ComponentListsCache();

  /** \brief Return the LocalDofLists object describing the DOFs corresponding
   * to
   *  AHMED matrix indices [start, start + indexCount). */
  shared_ptr<const ComponentLists> get(int start, int indexCount);

private:
  void findComponents(int start, int indexCount, std::vector<int> &pointIndices,
                      std::vector<std::vector<int>> &componentIndices,
                      std::vector<std::vector<int>> &arrayIndices) const;

  void findComponents(int index, std::vector<int> &pointIndices,
                      std::vector<std::vector<int>> &componentIndices,
                      std::vector<std::vector<int>> &arrayIndices) const;

private:
  /** \cond PRIVATE */
  const std::vector<std::size_t> &m_p2o;
  int m_componentCount;

  typedef tbb::concurrent_unordered_map<
      std::pair<int, int>, const ComponentLists *> ComponentListsMap;
  ComponentListsMap m_map;
  /** \endcond */
};

} // namespace Bempp

#endif
