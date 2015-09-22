#ifndef BEMPP_PY_SPHERE_HPP
#define BEMPP_PY_SPHERE_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <utility>
#include <cmath>

namespace Bempp {

class SphereMesh {
public:
  SphereMesh(int n, double radius = 1,
             std::array<double, 3> offset = {{0, 0, 0}});
  inline std::vector<std::array<double, 3>>* nodes() {
    return &m_nodes;
}

  inline std::vector<std::array<int, 3>>* elements() {
    return &m_elements;
    }

private:
  typedef std::map<std::pair<int, int>, int> LinesToVerticesMap;

  void refineTriangles();

  std::vector<std::array<double, 3>> m_nodes;
  std::vector<std::array<int, 3>> m_elements;
};

SphereMesh::SphereMesh(int n, double radius, std::array<double, 3> offset) {
  // Create initial octahedron

  for (int i = 0; i < 3; ++i) {
    std::array<double, 3> tmp{{0, 0, 0}};
    tmp[i] = 1;
    m_nodes.push_back(tmp);
    tmp[i] = -1;
    m_nodes.push_back(tmp);
  }

  m_elements.push_back(std::array<int, 3>({{2, 4, 0}}));
  m_elements.push_back(std::array<int, 3>({{1, 4, 2}}));
  m_elements.push_back(std::array<int, 3>({{3, 4, 1}}));
  m_elements.push_back(std::array<int, 3>({{0, 4, 3}}));
  m_elements.push_back(std::array<int, 3>({{5, 2, 0}}));
  m_elements.push_back(std::array<int, 3>({{5, 1, 2}}));
  m_elements.push_back(std::array<int, 3>({{5, 3, 1}}));
  m_elements.push_back(std::array<int, 3>({{5, 0, 3}}));

  // Now refine triangles as often as specified

  for (int i = 0; i < n; ++i)
    refineTriangles();

  // Apply offset and radius

  for (auto& node : m_nodes)
    for (int i = 0;i<3;++i) node[i] =radius*node[i]+offset[i];
}

void SphereMesh::refineTriangles() {

  int n = m_nodes.size();

  LinesToVerticesMap linesToVerticesMap;

  std::vector<std::array<int, 3>> newElements;

  for (const auto &elem : m_elements) {

    std::array<int, 3> newNodeNumbers;
    for (int i = 0; i < 3; ++i) {
      int p1 = i;
      int p2 = (i + 1) % 3;
      std::pair<int, int> t;
      if (elem[p1] < elem[p2])
        t = std::pair<int, int>(elem[p1], elem[p2]);
      else
        t = std::pair<int, int>(elem[p2], elem[p1]);

      auto search = linesToVerticesMap.find(t);
      if (search == end(linesToVerticesMap)) {
        linesToVerticesMap.insert(std::pair<std::pair<int, int>, int>(t, n));
        std::array<double, 3> newNode;
        newNode[0] = .5 * (m_nodes[elem[p1]][0] + m_nodes[elem[p2]][0]);
        newNode[1] = .5 * (m_nodes[elem[p1]][1] + m_nodes[elem[p2]][1]);
        newNode[2] = .5 * (m_nodes[elem[p1]][2] + m_nodes[elem[p2]][2]);
        newNodeNumbers[i] = n;
        m_nodes.push_back(newNode);
        ++n;
      } else
        newNodeNumbers[i] = search->second;
    }
    newElements.push_back(
        std::array<int, 3>({{elem[0], newNodeNumbers[0], newNodeNumbers[2]}}));
    newElements.push_back(
        std::array<int, 3>({{newNodeNumbers[0], elem[1], newNodeNumbers[1]}}));
    newElements.push_back(
        std::array<int, 3>({{newNodeNumbers[1], elem[2], newNodeNumbers[2]}}));
    newElements.push_back(std::array<int, 3>(
        {{newNodeNumbers[0], newNodeNumbers[1], newNodeNumbers[2]}}));
  }

  m_elements.clear();
  m_elements = newElements;

  for (auto &node : m_nodes) {
    double len = 0;
    for (int i = 0; i < 3; ++i)
      len += node[i] * node[i];
    len = std::sqrt(len);
    for (int i = 0; i < 3; ++i)
      node[i] = node[i] / len;
  }
}

}

#endif
