// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_DOF_PERMUTATION_IMPL_HPP
#define HMAT_DOF_PERMUTATION_IMPL_HPP

#include "dof_permutation.hpp"

namespace hmat {

inline DofPermutation::DofPermutation(std::size_t numberOfDofs)
    : m_numberOfDofs(numberOfDofs), m_originalDofToHMatDofMap(numberOfDofs),
      m_hMatDofToOriginalDofMap(numberOfDofs) {}

inline std::size_t
DofPermutation::mapOriginalDofToHMatDof(std::size_t originalDofIndex) const {

  return m_originalDofToHMatDofMap[originalDofIndex];
}

inline std::size_t
DofPermutation::mapHMatDofToOriginalDof(std::size_t hMatDofIndex) const {

  return m_hMatDofToOriginalDofMap[hMatDofIndex];
}

inline std::size_t DofPermutation::numberOfDofs() const {
  return m_numberOfDofs;
}

inline void DofPermutation::addDofIndexPair(std::size_t originalDofIndex,
                                            std::size_t hMatDofIndex) {

  m_originalDofToHMatDofMap[originalDofIndex] = hMatDofIndex;
  m_hMatDofToOriginalDofMap[hMatDofIndex] = originalDofIndex;
}

inline const std::vector<std::size_t> &
DofPermutation::hMatDofToOriginalDofMap() const {
  return m_hMatDofToOriginalDofMap;
}

inline const std::vector<std::size_t> &
DofPermutation::originalDofToHMatDofMap() const {
  return m_originalDofToHMatDofMap;
}
}

#endif
