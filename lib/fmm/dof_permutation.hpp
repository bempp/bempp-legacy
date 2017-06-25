#ifndef BEMPP_FMM_DOF_PERMUTATION_HPP
#define BEMPP_FMM_DOF_PERMUTATION_HPP

#include "fmm_common.hpp"
#include <vector>

namespace Fmm {

class DofPermutation {

public:
    DofPermutation(const std::vector<std::size_t>& p2o);

    std::size_t numberOfDofs() const;

    template <typename ValueType>
    void permute(const Eigen::Ref<const Vector<ValueType> >& vIn,
        Eigen::Ref<Vector<ValueType> >& vOut);

    template <typename ValueType>
    void unpermute(const Eigen::Ref<const Vector<ValueType> >& vIn,
        Eigen::Ref<Vector<ValueType> >& vOut);

private:
    std::vector<std::size_t> m_map;
};

} // namespace fmm

#include "dof_permutation_impl.hpp"

#endif
