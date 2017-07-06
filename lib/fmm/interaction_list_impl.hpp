#ifndef bempp_fmm_interaction_list_impl_hpp
#define bempp_fmm_interaction_list_impl_hpp

#include "interaction_list.hpp"
#include "unordered_set"

namespace Fmm {

inline InteractionList::InteractionList(const Octree &octree,
                                        unsigned long index, unsigned int level)
    : m_octree(octree), m_index(index), m_level(level), m_current(0) {

  if (m_level > 1) {

    std::vector<unsigned long> parentNeighbors;
    std::vector<unsigned long> indexNeighbors;
    // Get parent neighbors
    m_octree.getNeighbors(parentNeighbors, octree.getParent(m_index),
                          level - 1);
    // Get neighbors of index as a set
    m_octree.getNeighbors(indexNeighbors, m_index, level);
    std::unordered_set<unsigned long> neighborSet(begin(indexNeighbors),
                                                  end(indexNeighbors));

    // Collect all children of parent neighbors, which are not empty and
    // not a direct neighbor of the index.

    for (const auto &neighborIndex : parentNeighbors) {
      unsigned long firstChild = m_octree.getFirstChild(neighborIndex);
      unsigned long lastChild = m_octree.getLastChild(neighborIndex);
      for (unsigned long index = firstChild; index <= lastChild; index++)
        if (!m_octree.isEmpty(index, m_level) &&
            (neighborSet.count(index) == 0))
          m_interaction.push_back(index);
    }
  }
}

inline bool InteractionList::finished() const {
  return m_current >= m_interaction.size();
}

inline void InteractionList::next() { m_current++; }

inline const unsigned long &InteractionList::operator*() const {

  return m_interaction[m_current];
}
}

#endif
