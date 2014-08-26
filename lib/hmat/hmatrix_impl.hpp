// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_IMPL_HPP
#define HMAT_HMATRIX_IMPL_HPP

#include "hmatrix.hpp"
#include "hmatrix_data.hpp"

namespace hmat {

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree)
    : m_blockClusterTree(blockClusterTree) {}
}

#endif
