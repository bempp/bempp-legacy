// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_HPP
#define HMAT_HMATRIX_HPP

#include "common.hpp"
#include "block_cluster_tree.hpp"
#include "hmatrix_compressor.hpp"
#include "data_accessor.hpp"
#include <armadillo>
#include <unordered_map>

namespace hmat {

template <typename ValueType> class HMatrixData;
template <typename ValueType, int N> class HMatrix;

template <typename ValueType>
using DefaultHMatrixType = HMatrix<ValueType,2>;

template <typename ValueType, int N> class HMatrix {
public:
  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree);
  HMatrix(const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
          const HMatrixCompressor<ValueType,N> &hMatrixCompressor,
          const DataAccessor<ValueType,N> &dataAccessor);
  
  void initialize(const HMatrixCompressor<ValueType,N> &hMatrixCompressor,
                  const DataAccessor<ValueType,N> &dataAccessor);

  void reset();

private:
  shared_ptr<BlockClusterTree<N>> m_blockClusterTree;
  std::unordered_map<shared_ptr<BlockClusterTreeNode<N>>,
                     shared_ptr<HMatrixData<ValueType>>> m_hMatrixData;
};

}

#include "hmatrix_impl.hpp"

#endif // HMATRIX_HPP
