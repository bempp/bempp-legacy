#ifndef bempp_fmm_chebychev_tools_impl_hpp
#define bempp_fmm_chebychev_tools_impl_hpp

#include "chebychev_tools.hpp"

#include <cmath>

namespace Fmm {

inline ChebychevTools::ChebychevTools(int order) : m_terms(1 + order) {

  m_nodes.resize(m_terms);
  for (int m = 0; m < m_terms; ++m)
    m_nodes(m) = cos(M_PI * double(2 * m + 1) / double(2 * m_terms));

  m_chebychevPolValuesAtNodes.resize(m_terms, m_terms);
  // compute T_k(x) for k between 0 and N-1 inclusive.
  for (unsigned int m = 0; m < m_terms; m++) {
    m_chebychevPolValuesAtNodes(m, 0) = 1;
    m_chebychevPolValuesAtNodes(m, 1) = m_nodes[m];
    for (unsigned int k = 2; k < m_terms; k++)
      m_chebychevPolValuesAtNodes(m, k) =
          2 * m_nodes[m] * m_chebychevPolValuesAtNodes(m, k - 1) -
          m_chebychevPolValuesAtNodes(m, k - 2);
  }
}

inline const Vector<double> &ChebychevTools::chebychevNodes() const {
  return m_nodes;
}

const Matrix<double> &ChebychevTools::chebychevPolValuesAtNodes() const {

  return m_chebychevPolValuesAtNodes;
}
}

#endif
