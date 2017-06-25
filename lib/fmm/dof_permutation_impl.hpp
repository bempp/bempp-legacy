#ifndef BEMPP_FMM_DOF_PERMUTATION_IMPL_HPP
#define BEMPP_FMM_DOF_PERMUTATION_IMPL_HPP

#include "dof_permutation.hpp"

namespace Fmm {

inline DofPermutation::DofPermutation(const std::vector<std::size_t>& p2o)
    : m_map(p2o)
{
}

inline std::size_t DofPermutation::numberOfDofs() const { return m_map.size(); }

template <typename ValueType>
inline void
DofPermutation::unpermute(const Eigen::Ref<const Vector<ValueType> >& vIn,
    Eigen::Ref<Vector<ValueType> >& vOut)
{
    for (std::size_t i = 0; i < numberOfDofs(); ++i)
        vOut(m_map[i]) = vIn(i);
}

template <typename ValueType>
inline void
DofPermutation::permute(const Eigen::Ref<const Vector<ValueType> >& vIn,
    Eigen::Ref<Vector<ValueType> >& vOut)
{
    for (std::size_t i = 0; i < numberOfDofs(); ++i)
        vOut(i) = vIn(m_map[i]);
}

} // namespace fmm

#endif
