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

#include "local_dof_lists_cache.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"

#include <map>
#include <set>
#include <utility>

namespace Bempp {

template <typename BasisFunctionType>
LocalDofListsCache<BasisFunctionType>::LocalDofListsCache(
    const Space<BasisFunctionType> &space, const std::vector<std::size_t> &p2o,
    bool indexWithGlobalDofs)
    : m_space(space), m_p2o(p2o), m_indexWithGlobalDofs(indexWithGlobalDofs) {}

template <typename BasisFunctionType>
LocalDofListsCache<BasisFunctionType>::~LocalDofListsCache() {
  for (typename LocalDofListsMap::const_iterator it = m_map.begin();
       it != m_map.end(); ++it)
    delete it->second;
  m_map.clear();
}

template <typename BasisFunctionType>
shared_ptr<const LocalDofLists<BasisFunctionType>>
LocalDofListsCache<BasisFunctionType>::get(int start, int indexCount) {
  if (indexCount == 1) {
    shared_ptr<LocalDofLists<BasisFunctionType>> result(
        new LocalDofLists<BasisFunctionType>);
    findLocalDofs(start, result->originalIndices, result->elementIndices,
                  result->localDofIndices, result->localDofWeights,
                  result->arrayIndices);
    return result;
  }

  std::pair<int, int> key(start, indexCount);
  typename LocalDofListsMap::const_iterator it = m_map.find(key);
  if (it != m_map.end()) {
    return make_shared_from_ref(*it->second);
  }

  // The relevant local DOF list doesn't exist yet and must be created.
  LocalDofLists<BasisFunctionType> *newLists =
      new LocalDofLists<BasisFunctionType>;
  findLocalDofs(start, indexCount, newLists->originalIndices,
                newLists->elementIndices, newLists->localDofIndices,
                newLists->localDofWeights, newLists->arrayIndices);

  // Attempt to insert the newly created DOF list into the map
  std::pair<typename LocalDofListsMap::iterator, bool> result =
      m_map.insert(std::make_pair(key, newLists));
  if (result.second)
    // Insertion succeeded. The newly created DOF list will be deleted in
    // our own destructor
    ;
  else
    // Insertion failed -- another thread was faster. Delete the newly
    // created DOF list.
    delete newLists;

  // Return pointer to the DOF list that ended up in the map.
  return make_shared_from_ref(*result.first->second);
}

template <typename BasisFunctionType>
void LocalDofListsCache<BasisFunctionType>::findLocalDofs(
    int start, int indexCount,
    std::vector<typename LocalDofLists<BasisFunctionType>::DofIndex> &
        originalIndices,
    std::vector<int> &elementIndices,
    std::vector<std::vector<LocalDofIndex>> &localDofIndices,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights,
    std::vector<std::vector<int>> &arrayIndices) const {
  using std::make_pair;
  using std::map;
  using std::pair;
  using std::set;
  using std::vector;

  // Convert permuted indices into original indices
  originalIndices.resize(indexCount);
  for (int i = 0; i < indexCount; ++i)
    originalIndices[i] = m_p2o[start + i];

  // map of pairs (local dof index, array index) to local dof weights
  typedef std::map<pair<LocalDofIndex, int>, BasisFunctionType>
      LocalDofWeightMap;
  typedef std::map<EntityIndex, LocalDofWeightMap> LocalDofMap;

  // Temporary map: entityIndex -> set(localDofIndex, arrayIndex)
  // with arrayIndex standing for the index of the row or column in the matrix
  // that needs to be returned to Ahmed.
  LocalDofMap requiredLocalDofs;

  // Retrieve lists of local DOFs corresponding to original indices,
  // treated either as global DOFs (if m_indexWithGlobalDofs is true)
  // or flat local DOFs (if m_indexWithGlobalDofs is false)
  if (m_indexWithGlobalDofs) {
    vector<vector<LocalDof>> localDofs;
    vector<vector<BasisFunctionType>> localDofWeights;
    m_space.global2localDofs(originalIndices, localDofs, localDofWeights);

    for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex) {
      const vector<LocalDof> &currentLocalDofs = localDofs[arrayIndex];
      const vector<BasisFunctionType> &currentLocalDofWeights =
          localDofWeights[arrayIndex];
      for (size_t j = 0; j < currentLocalDofs.size(); ++j)
        requiredLocalDofs[currentLocalDofs[j].entityIndex][make_pair(
            currentLocalDofs[j].dofIndex, arrayIndex)] =
            currentLocalDofWeights[j];
    }
  } else {
    vector<LocalDof> localDofs;
    m_space.flatLocal2localDofs(originalIndices, localDofs);

    for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex) {
      const LocalDof &currentLocalDof = localDofs[arrayIndex];
      requiredLocalDofs[currentLocalDof.entityIndex][make_pair(
          currentLocalDof.dofIndex, arrayIndex)] = 1.;
    }
  }

  // Use the temporary map requiredLocalDofs to build the three output vectors
  const int elementCount = requiredLocalDofs.size();
  // vector<EntityIndex> elementIndices;
  // elementIndices.reserve(elementCount);

  elementIndices.resize(elementCount);
  localDofIndices.clear();
  localDofIndices.resize(elementCount);
  localDofWeights.clear();
  localDofWeights.resize(elementCount);
  arrayIndices.clear();
  arrayIndices.resize(elementCount);

  // const ReverseElementMapper& mapper = m_view.reverseElementMapper();

  int e = 0;
  for (typename LocalDofMap::const_iterator mapIt = requiredLocalDofs.begin();
       mapIt != requiredLocalDofs.end(); ++mapIt, ++e) {
    elementIndices[e] = mapIt->first;
    //        elements[e] = &mapper.entityPointer(mapIt->first);
    for (typename LocalDofWeightMap::const_iterator wmapIt =
             mapIt->second.begin();
         wmapIt != mapIt->second.end(); ++wmapIt) {
      localDofIndices[e].push_back(wmapIt->first.first);
      localDofWeights[e].push_back(wmapIt->second);
      arrayIndices[e].push_back(wmapIt->first.second);
    }
  }
}

template <typename BasisFunctionType>
void LocalDofListsCache<BasisFunctionType>::findLocalDofs(
    int index,
    std::vector<typename LocalDofLists<BasisFunctionType>::DofIndex> &
        originalIndices,
    std::vector<int> &elementIndices,
    std::vector<std::vector<LocalDofIndex>> &localDofIndices,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights,
    std::vector<std::vector<int>> &arrayIndices) const {
  using std::make_pair;
  using std::map;
  using std::pair;
  using std::set;
  using std::vector;

  // Convert permuted indices into original indices
  originalIndices.resize(1);
  assert(index >= 0 && index < m_p2o.size());
  originalIndices[0] = m_p2o[index];

  // set of pairs (local dof index, array index)
  typedef std::set<pair<LocalDofIndex, int>> LocalDofSet;
  typedef std::map<EntityIndex, LocalDofSet> LocalDofMap;

  // Retrieve lists of local DOFs corresponding to original indices,
  // treated either as global DOFs (if m_indexWithGlobalDofs is true)
  // or flat local DOFs (if m_indexWithGlobalDofs is false)
  if (m_indexWithGlobalDofs) {
    // raw -- means arrays as returned by global2localDofs; the data need to be
    // subsequently rearranged
    vector<vector<LocalDof>> rawLocalDofs;
    vector<vector<BasisFunctionType>> rawLocalDofWeights;
    m_space.global2localDofs(originalIndices, rawLocalDofs, rawLocalDofWeights);

    // Here we assume that no global DOF contains more than one local DOF
    // from a particular element
    const vector<LocalDof> &currentLocalDofs = rawLocalDofs[0];
    const vector<BasisFunctionType> &currentLocalDofWeights =
        rawLocalDofWeights[0];
    size_t cnt = currentLocalDofs.size();
    elementIndices.resize(cnt);
    localDofIndices.resize(cnt, std::vector<LocalDofIndex>(1));
    localDofWeights.resize(cnt, std::vector<BasisFunctionType>(1));
    arrayIndices.resize(cnt, std::vector<int>(1));
    for (size_t j = 0; j < cnt; ++j) {
      elementIndices[j] = currentLocalDofs[j].entityIndex;
      localDofIndices[j][0] = currentLocalDofs[j].dofIndex;
      localDofWeights[j][0] = currentLocalDofWeights[j];
      arrayIndices[j][0] = 0;
    }
  } else {
    vector<LocalDof> localDofs;
    m_space.flatLocal2localDofs(originalIndices, localDofs);

    size_t cnt = localDofs.size();
    elementIndices.resize(cnt);
    localDofIndices.resize(cnt, std::vector<LocalDofIndex>(1));
    localDofWeights.resize(cnt, std::vector<BasisFunctionType>(1));
    arrayIndices.resize(cnt, std::vector<int>(1));

    for (size_t j = 0; j < cnt; ++j) {
      elementIndices[j] = localDofs[0].entityIndex;
      localDofIndices[j][0] = localDofs[0].dofIndex;
      localDofWeights[j][0] = 1.;
      arrayIndices[j][0] = 0;
    }
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(LocalDofListsCache);

} // namespace Bempp
