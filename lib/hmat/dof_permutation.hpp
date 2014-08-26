// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_DOF_PERMUTATION_HPP
#define HMAT_DOF_PERMUTATION_HPP

#include "common.hpp"
#include <vector>

namespace hmat {

class DofPermutation {

public:
  DofPermutation(std::size_t numberOfDofs);

  std::size_t mapOriginalDofToHMatDof(std::size_t originalDofIndex) const;
  std::size_t mapHMatDofToOriginalDof(std::size_t hMatDofIndex) const;
  std::size_t numberOfDofs() const;

  void addDofIndexPair(std::size_t originalDofIndex, std::size_t hMatDofIndex);

private:
  std::size_t m_numberOfDofs;
  std::vector<std::size_t> m_originalDofToHMatDofMap;
  std::vector<std::size_t> m_hMatDofToOriginalDofMap;
};
}

#include "dof_permutation_impl.hpp"

#endif
