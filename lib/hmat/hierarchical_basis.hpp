
#ifndef HMAT_HIERARCHICAL_BASIS_HPP
#define HMAT_HIERARCHICAL_BASIS_HPP

#include "common.hpp"
#include "cluster_tree.hpp"
#include "tbb/concurrent_unordered_map.h"


template <typename ValueType> class TransferOperator;

template <typename ValueType, int N>
class HierarchicalBasis:
{

public:

  typedef tbb::concurrent_unordered_map<
      shared_ptr<const ClusterTreeNode<N>>, shared_ptr<TransferOperator<ValueType>>,
      shared_ptr_hash<ClusterTreeNode<N>>> TransferOperatorContainer;

  typedef tbb::concurrent_unordered_map<
      shared_ptr<const ClusterTreeNode<N>>, shared_ptr<LeafBasis<ValueType>>,
      shared_ptr_hash<ClusterTreeNode<N>>> LeafBasisContainer;

    HierarchicalBasis(const shared_ptr<const ClusterTree<N>>& clusterTree);
    virtual ~HierarchicalBasis();

    void addTransferOperator(
            const shared_ptr<const ClusterTreeNode<N>& node,
            const shared_ptr<const TransferOperator>& transferOperator);
    void addLeafBasis(
            const shared_ptr<const ClusterTreeNode<N>& node,
            const shared_ptr<const LeafBasis>& leafBasis);


private:

    shared_ptr<const ClusterTree<N>>& m_clusterTree;
    
    TransferOperatorContainer m_transferOperatorContainer;
    LeafBasisContainer m_leafBasisContainer;


};

#include "hierarchical_basis_impl.hpp"


#endif
