// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_SIMPLE_TREE_NODE_HPP
#define HMAT_SIMPLE_TREE_NODE_HPP

#include "common.hpp"
#include <array>
#include <memory>

namespace hmat {
template <typename T, int N>
class SimpleTreeNode : public enable_shared_from_this<SimpleTreeNode<T, N>> {
public:
  SimpleTreeNode(const T &data);
  SimpleTreeNode(const shared_ptr<SimpleTreeNode<T, N>> &root, const T &data);
  const shared_ptr<const SimpleTreeNode<T, N>> root() const;
  const shared_ptr<const SimpleTreeNode<T, N>> child(int i) const;
  const shared_ptr<SimpleTreeNode<T, N>> child(int i);

  const T &data() const;
  T &data();

  void addChild(const T &child, int i);
  void removeChildren();

  void addSubTree(shared_ptr<SimpleTreeNode<T, N>> &subTree, int i);

  bool isLeaf() const;

  const std::vector<shared_ptr<const SimpleTreeNode<T, N>>> leafNodes() const;
  const std::vector<shared_ptr<SimpleTreeNode<T, N>>> leafNodes();

  std::size_t numberOfLeafs() const;

private:
  std::array<shared_ptr<SimpleTreeNode<T, N>>, N> m_children;
  weak_ptr<SimpleTreeNode<T, N>> m_root;

  T m_data;
};
}

#include "simple_tree_node_impl.hpp"

#endif
