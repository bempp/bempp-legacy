// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_SIMPLE_TREE_NODE_IMPL_HPP
#define HMAT_SIMPLE_TREE_NODE_IMPL_HPP

#include "common.hpp"
#include <array>
#include <stdexcept>
#include <cassert>
#include <functional>

#include "simple_tree_node.hpp"

namespace hmat {
template <typename T, int N>
SimpleTreeNode<T, N>::SimpleTreeNode(const T &data)
    : m_root(shared_ptr<SimpleTreeNode<T, N>>()), m_data(data){};

template <typename T, int N>
SimpleTreeNode<T, N>::SimpleTreeNode(
    const shared_ptr<SimpleTreeNode<T, N>> &root, const T &data)
    : m_root(root), m_data(data){};

template <typename T, int N>
const shared_ptr<const SimpleTreeNode<T, N>>
SimpleTreeNode<T, N>::root() const {
  return m_root.lock();
}

template <typename T, int N>
const shared_ptr<const SimpleTreeNode<T, N>>
SimpleTreeNode<T, N>::child(int i) const {
  assert(i < N);
  assert(m_children[i]);

  return m_children[i];
}

template <typename T, int N>
const shared_ptr<SimpleTreeNode<T, N>> SimpleTreeNode<T, N>::child(int i) {
  assert(i < N);
  assert(m_children[i]);

  return m_children[i];
}

template <typename T, int N> const T &SimpleTreeNode<T, N>::data() const {
  return m_data;
}

template <typename T, int N> T &SimpleTreeNode<T, N>::data() { return m_data; }

template <typename T, int N>
void SimpleTreeNode<T, N>::addChild(const T &data, int i) {

  assert(i < N);

  m_children[i] =
      make_shared<SimpleTreeNode<T, N>>(this->shared_from_this(), data);
}

template <typename T, int N>
void SimpleTreeNode<T, N>::removeChildren() {

  for (auto& child: m_children)
    child.reset();

}

template <typename T, int N>
void SimpleTreeNode<T, N>::addSubTree(shared_ptr<SimpleTreeNode<T, N>> &subTree,
                                      int i) {

  subTree->m_root = this->shared_from_this();
  m_children[i] = subTree;
}

template <typename T, int N> bool SimpleTreeNode<T, N>::isLeaf() const {

  for (auto child : m_children)
    if (child)
      return false;

  return true;
}

template <typename T, int N>
const std::vector<shared_ptr<const SimpleTreeNode<T, N>>>
SimpleTreeNode<T, N>::leafNodes() const {

  std::function<void(SimpleTreeNode<T, N> &)> getLeafsImpl;

  std::vector<shared_ptr<const SimpleTreeNode<T, N>>> leafVector;

  getLeafsImpl = [&leafVector, &getLeafsImpl](SimpleTreeNode<T, N> &node) {

    if (node.isLeaf())
      leafVector.push_back(node);
    else
      for (int i = 0; i < N; ++i)
        getLeafsImpl(*(node.child(i)));
  };

  getLeafsImpl(*this);
  return leafVector;
}

template <typename T, int N>
const std::vector<shared_ptr<SimpleTreeNode<T, N>>>
SimpleTreeNode<T, N>::leafNodes() {

  std::function<void(SimpleTreeNode<T, N> &)> getLeafsImpl;

  std::vector<shared_ptr<SimpleTreeNode<T, N>>> leafVector;

  getLeafsImpl = [&leafVector, &getLeafsImpl](SimpleTreeNode<T, N> &node) {

    if (node.isLeaf())
      leafVector.push_back(node.shared_from_this());
    else
      for (int i = 0; i < N; ++i)
        getLeafsImpl(*(node.child(i)));
  };

  getLeafsImpl(*this);
  return leafVector;
}
}
#endif
