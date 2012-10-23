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

#include <auto_ptr.h>
#include <map>
#include <set>
#include <utility>

namespace Bempp
{

template <typename BasisFunctionType>
LocalDofListsCache<BasisFunctionType>::LocalDofListsCache(
        const Space<BasisFunctionType>& space,
        const std::vector<unsigned int>& p2o,
        bool indexWithGlobalDofs) :
    m_space(space), m_p2o(p2o), m_indexWithGlobalDofs(indexWithGlobalDofs)
{
}

template <typename BasisFunctionType>
shared_ptr<const LocalDofLists> LocalDofListsCache<BasisFunctionType>::get(
        int start, int indexCount)
{
    if (indexCount == 1) {
        shared_ptr<LocalDofLists> result(new LocalDofLists);
        findLocalDofs(start, result->originalIndices, result->elementIndices,
                      result->localDofIndices, result->arrayIndices);
        return result;
    }

    std::pair<int, int> key(start, indexCount);
    typename LocalDofListsMap::const_iterator it = m_map.find(key);
    if (it != m_map.end()) {
        return it->second;
    }

    // The relevant local DOF list doesn't exist yet and must be created.
    shared_ptr<LocalDofLists> newLists(new LocalDofLists);
    findLocalDofs(start, indexCount,
                  newLists->originalIndices, newLists->elementIndices,
                  newLists->localDofIndices, newLists->arrayIndices);

    // Attempt to insert the newly created DOF list into the map
    std::pair<typename LocalDofListsMap::iterator, bool> result =
            m_map.insert(std::make_pair(key, newLists));

    // Return pointer to the DOF list that ended up in the map.
    return result.first->second;
}

template <typename BasisFunctionType>
void LocalDofListsCache<BasisFunctionType>::
findLocalDofs(
        int start,
        int indexCount,
        std::vector<LocalDofLists::DofIndex>& originalIndices,
        std::vector<int>& elementIndices,
        std::vector<std::vector<LocalDofIndex> >& localDofIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    using std::make_pair;
    using std::map;
    using std::pair;
    using std::set;
    using std::vector;

    // Convert permuted indices into original indices
    originalIndices.resize(indexCount);
    for (int i = 0; i < indexCount; ++i)
        originalIndices[i] = m_p2o[start + i];

    // set of pairs (local dof index, array index)
    typedef std::set<pair<LocalDofIndex, int> > LocalDofSet;
    typedef std::map<EntityIndex, LocalDofSet> LocalDofMap;

    // Temporary map: entityIndex -> set(localDofIndex, arrayIndex)
    // with arrayIndex standing for the index of the row or column in the matrix
    // that needs to be returned to Ahmed.
    LocalDofMap requiredLocalDofs;

    // Retrieve lists of local DOFs corresponding to original indices,
    // treated either as global DOFs (if m_indexWithGlobalDofs is true)
    // or flat local DOFs (if m_indexWithGlobalDofs is false)
    if (m_indexWithGlobalDofs)
    {
        vector<vector<LocalDof> > localDofs;
        m_space.global2localDofs(originalIndices, localDofs);

        for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex)
        {
            const vector<LocalDof>& currentLocalDofs = localDofs[arrayIndex];
            for (size_t j = 0; j < currentLocalDofs.size(); ++j)
                requiredLocalDofs[currentLocalDofs[j].entityIndex]
                        .insert(make_pair(currentLocalDofs[j].dofIndex, arrayIndex));
        }
    }
    else
    {
        vector<LocalDof> localDofs;
        m_space.flatLocal2localDofs(originalIndices, localDofs);

        for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex)
        {
            const LocalDof& currentLocalDof = localDofs[arrayIndex];
            requiredLocalDofs[currentLocalDof.entityIndex]
                    .insert(make_pair(currentLocalDof.dofIndex, arrayIndex));
        }
    }

    // Use the temporary map requiredLocalDofs to build the three output vectors
    const int elementCount = requiredLocalDofs.size();
    // vector<EntityIndex> elementIndices;
    // elementIndices.reserve(elementCount);

    elementIndices.resize(elementCount);
    localDofIndices.clear();
    localDofIndices.resize(elementCount);
    arrayIndices.clear();
    arrayIndices.resize(elementCount);

    // const ReverseElementMapper& mapper = m_view.reverseElementMapper();

    int e = 0;
    for (LocalDofMap::const_iterator mapIt = requiredLocalDofs.begin();
         mapIt != requiredLocalDofs.end(); ++mapIt, ++e)
    {
        elementIndices[e] = mapIt->first;
//        elements[e] = &mapper.entityPointer(mapIt->first);
        for (LocalDofSet::const_iterator setIt = mapIt->second.begin();
             setIt != mapIt->second.end(); ++setIt)
        {
            localDofIndices[e].push_back(setIt->first);
            arrayIndices[e].push_back(setIt->second);
        }
    }
}

template <typename BasisFunctionType>
void LocalDofListsCache<BasisFunctionType>::
findLocalDofs(
        int index,
        std::vector<LocalDofLists::DofIndex>& originalIndices,
        std::vector<int>& elementIndices,
        std::vector<std::vector<LocalDofIndex> >& localDofIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    using std::make_pair;
    using std::map;
    using std::pair;
    using std::set;
    using std::vector;

    // Convert permuted indices into original indices
    originalIndices.resize(1);
    originalIndices[0] = m_p2o[index];

    // set of pairs (local dof index, array index)
    typedef std::set<pair<LocalDofIndex, int> > LocalDofSet;
    typedef std::map<EntityIndex, LocalDofSet> LocalDofMap;

    // Retrieve lists of local DOFs corresponding to original indices,
    // treated either as global DOFs (if m_indexWithGlobalDofs is true)
    // or flat local DOFs (if m_indexWithGlobalDofs is false)
    if (m_indexWithGlobalDofs) {
        vector<vector<LocalDof> > localDofs;
        m_space.global2localDofs(originalIndices, localDofs);

        // Here we assume that no global DOF contains more than one local DOF
        // from a particular element
        const vector<LocalDof>& currentLocalDofs = localDofs[0];
        size_t cnt = currentLocalDofs.size();
        elementIndices.resize(cnt);
        localDofIndices.resize(cnt, std::vector<LocalDofIndex>(1));
        arrayIndices.resize(cnt, std::vector<int>(1));
        for (size_t j = 0; j < cnt; ++j) {
            elementIndices[j] = currentLocalDofs[j].entityIndex;
            localDofIndices[j][0] = currentLocalDofs[j].dofIndex;
            arrayIndices[j][0] = 0;
        }
    } else {
        vector<LocalDof> localDofs;
        m_space.flatLocal2localDofs(originalIndices, localDofs);

        size_t cnt = localDofs.size();
        elementIndices.resize(cnt);
        localDofIndices.resize(cnt, std::vector<LocalDofIndex>(1));
        arrayIndices.resize(cnt, std::vector<int>(1));

        for (size_t j = 0; j < cnt; ++j) {
            elementIndices[j] = localDofs[0].entityIndex;
            localDofIndices[j][0] = localDofs[0].dofIndex;
            arrayIndices[j][0] = 0;
        }
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(LocalDofListsCache);

} // namespace Bempp
