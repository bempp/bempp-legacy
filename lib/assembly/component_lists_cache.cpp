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

#include "component_lists_cache.hpp"

#include <map>
#include <set>
#include <utility>

namespace Bempp
{

ComponentListsCache::ComponentListsCache(
        const std::vector<unsigned int>& p2o, int componentCount) :
    m_p2o(p2o), m_componentCount(componentCount)
{
}

ComponentListsCache::~ComponentListsCache()
{
    for (ComponentListsMap::const_iterator it = m_map.begin();
         it != m_map.end(); ++it)
        delete it->second;
    m_map.clear();
}

shared_ptr<const ComponentLists>
ComponentListsCache::get(int start, int indexCount)
{
    if (indexCount == 1) {
        shared_ptr<ComponentLists > result(
            new ComponentLists);
        findComponents(start, result->pointIndices,
                       result->componentIndices, result->arrayIndices);
        return result;
    }

    std::pair<int, int> key(start, indexCount);
    ComponentListsMap::const_iterator it = m_map.find(key);
    if (it != m_map.end()) {
        return make_shared_from_ref(*it->second);
    }

    // The relevant local DOF list doesn't exist yet and must be created.
    ComponentLists* newLists = new ComponentLists;
    findComponents(start, indexCount, newLists->pointIndices,
                   newLists->componentIndices, newLists->arrayIndices);

    // Attempt to insert the newly created DOF list into the map
    std::pair<ComponentListsMap::iterator, bool> result =
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

void ComponentListsCache::findComponents(
        int start,
        int indexCount,
        std::vector<int>& pointIndices,
        std::vector<std::vector<int> >& componentIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    // map of pairs (local dof index, array index) to local dof weights
    typedef int PointIndex;
    typedef int ComponentIndex;
    typedef int ArrayIndex;
    typedef std::pair<ComponentIndex, ArrayIndex> ComponentArrayIndexPair;
    typedef std::set<ComponentArrayIndexPair> ComponentSet;
    typedef std::map<PointIndex, ComponentSet> ComponentMap;

    // Temporary map: pointIndex -> set(componentIndex, arrayIndex)
    // with arrayIndex standing for the index of the row or column in the matrix
    // that needs to be returned to Ahmed.
    ComponentMap requiredComponents;

    // Retrieve lists of local DOFs corresponding to original indices,
    // treated either as global DOFs (if m_indexWithGlobalDofs is true)
    // or flat local DOFs (if m_indexWithGlobalDofs is false)
    for (int i = 0; i < indexCount; ++i) {
        PointIndex pointIndex = m_p2o[start + i] / m_componentCount;
        ComponentIndex componentIndex = m_p2o[start + i] % m_componentCount;
        requiredComponents[pointIndex].insert(std::make_pair(componentIndex, i));
    }

    // Use the temporary map requiredComponents to build the three output vectors
    const int pointCount = requiredComponents.size();
    // vector<EntityIndex> elementIndices;
    // elementIndices.reserve(elementCount);

    pointIndices.resize(pointCount);
    componentIndices.clear();
    componentIndices.resize(pointCount);
    arrayIndices.clear();
    arrayIndices.resize(pointCount);

    int p = 0;
    for (ComponentMap::const_iterator mapIt = requiredComponents.begin();
         mapIt != requiredComponents.end(); ++mapIt, ++p) {
        pointIndices[p] = mapIt->first;
        for (ComponentSet::const_iterator setIt = mapIt->second.begin();
             setIt != mapIt->second.end(); ++setIt) {
            componentIndices[p].push_back(setIt->first);
            arrayIndices[p].push_back(setIt->second);
        }
    }
}

void ComponentListsCache::findComponents(
        int index,
        std::vector<int>& pointIndices,
        std::vector<std::vector<int> >& componentIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    typedef int PointIndex;
    typedef int ComponentIndex;
    PointIndex pointIndex = m_p2o[index] / m_componentCount;
    ComponentIndex componentIndex = m_p2o[index] % m_componentCount;
    pointIndices.resize(1);
    pointIndices[0] = pointIndex;
    componentIndices.resize(1);
    componentIndices[0].resize(1);
    componentIndices[0][0] = componentIndex;
    arrayIndices.resize(1);
    arrayIndices[0].resize(1);
    arrayIndices[0][0] = 0;
}

} // namespace Bempp
