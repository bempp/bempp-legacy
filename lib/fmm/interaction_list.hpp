#ifndef bempp_fmm_interaction_list_hpp
#define bempp_fmm_interaction_list_hpp

#include "fmm_common.hpp"
#include "octree.hpp"
#include <vector>

namespace Fmm {

class InteractionList {

public:
  InteractionList(const Octree &octree, unsigned long index,
                  unsigned int level);

  /** \brief Return the next nonempty cube in the interaction list */
  void next();

  /** True if finished iterating */
  bool finished() const;

  /** Get reference to current element */
  const unsigned long &operator*() const;

private:
  std::vector<unsigned long> m_interaction;
  const Octree &m_octree;
  unsigned long m_index;
  unsigned int m_level;
  unsigned int m_current;
};
}

#include "interaction_list_impl.hpp"

#endif
