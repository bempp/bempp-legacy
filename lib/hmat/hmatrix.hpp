// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_HPP
#define HMAT_HMATRIX_HPP

#include "common.hpp"
#include "block_cluster_tree.hpp"
#include <armadillo>
#include <unordered_map>

namespace hmat {

template <typename ValueType> class HMatrixData;

template <typename ValueType, int N = 2> class HMatrix {

public:
  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree);

private:
  shared_ptr<BlockClusterTree<N>> m_blockClusterTree;
  std::unordered_map<shared_ptr<BlockClusterTreeNode<N>>,
                     shared_ptr<HMatrixData<ValueType>>> m_hMatrixData;
};
}

#include "hmatrix_impl.hpp"

#endif // HMATRIX_HPP
