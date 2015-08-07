#ifndef HMAT_HIERARCHICAL_BASIS_IMPL_HPP
#define HMAT_HIERARCHICAL_BASIS_IMPL_HPP

#include "hierarchical_basis.hpp"

template <typename ValueType, int N>
HierarchicalBasis<ValueType, N>::HierarchicalBasis(
        const shared_ptr<const ClusterTree<N>>& clusterTree):
    m_clusterTree(clusterTree){}

template <typename ValueType, int N>
HierarchicalBasis<ValueType, N>::~HierarchicalBasis(){}

template <typename ValueType, int N>
HierarchicalBasis<ValueType, N>::addTransferOperator(
        const shared_ptr<const ClusterTreeNode<N>> node,
        const shared_ptr<const TransferOperator>& transferOperator)
{
    m_transferOperatorContainer[node] = transferOperator;
}


template <typename ValueType, int N>
HierarchicalBasis<ValueType, N>::addLeafBasis(
        const shared_ptr<const ClusterTreeNode<N>> node,
        const shared_ptr<const TransferOperator>& leafBasis)
{
    if (!node->isLeaf())
        throw std::runtime_error(
                "HierarchicalBasis::addLeafBasis: "
                 "Node must be a leaf node.");
    m_leafBasisContainer[node] = leafBasis;
}



#endif
