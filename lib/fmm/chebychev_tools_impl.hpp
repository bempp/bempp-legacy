#ifndef bempp_fmm_chebychev_tools_impl_hpp
#define bempp_fmm_chebychev_tools_impl_hpp

#include "chebychev_tools.hpp"

#include <cmath>

namespace Fmm {

inline ChebychevTools::ChebychevTools(int order) : m_terms(1 + order) {

  assert(order > 0);
  m_nodes.resize(m_terms);
  for (int m = 0; m < m_terms; ++m)
    m_nodes(m) = cos(M_PI * double(m) / double(order));

  m_chebychevPolValuesAtNodes.resize(m_terms, m_terms);
  // compute T_k(x) for k between 0 and N-1 inclusive.

  m_chebychevPolValuesAtNodes.col(0).array() = 1;

  if (m_terms > 1) {
    m_chebychevPolValuesAtNodes.col(1).array() = m_nodes.array();

    for (int k = 2; k < m_terms; k++)
      m_chebychevPolValuesAtNodes.col(k).array() =
          2 * m_nodes.array() * m_chebychevPolValuesAtNodes.col(k - 1).array() -
          m_chebychevPolValuesAtNodes.col(k - 2).array();
  }

  m_barycentricWeights.resize(m_terms);
  int sign = 1;
  for (int k = 0; k < m_terms; ++k) {
    m_barycentricWeights[k] = sign;
    sign *= -1;
  }
  m_barycentricWeights[0] = double(1) / 2;
  m_barycentricWeights[order] *= double(1) / 2;
}

inline const Vector<double> &ChebychevTools::chebychevNodes() const {
  return m_nodes;
}

inline const Matrix<double> &ChebychevTools::chebychevPolValuesAtNodes() const {

  return m_chebychevPolValuesAtNodes;
}

inline void ChebychevTools::evaluateInterpolationPolynomial(
    const Vector<double> &weights, const Vector<double> &evaluationPoints,
    Vector<double> &result) const {

  int n = evaluationPoints.size();
  Vector<double> numer = Vector<double>::Zero(n);
  Vector<double> denom = Vector<double>::Zero(n);
  Vector<int> exact = Vector<int>::Zero(n);
  Vector<double> xdiff = Vector<double>::Zero(n);

  Vector<int> ones = Vector<int>::Ones(n);
  Vector<int> zeros = Vector<int>::Zero(n);

  for (int j = 0; j < m_terms; ++j) {
    xdiff.array() = evaluationPoints.array() - m_nodes(j);
    Vector<double> tmp = m_barycentricWeights(j) / xdiff.array();
    numer.array() += tmp.array() * weights(j);
    denom.array() += tmp.array();
    exact.array() =
        (xdiff.array() == 0).select((1 + j) * ones.array(), exact.array());
  }

  result = numer.array() / denom.array();

  for (int j = 0; j < n; ++j)
    if (exact(j) > 0)
      result(j) = weights(exact(j) - 1);
}

inline Matrix<double>
ChebychevTools::interpolateToChildren(double parentLength,
                                      double childLength) const {

  double r = 2 * childLength / parentLength;
  Vector<double> points(2 * m_terms);
  for (int j = 0; j < m_terms; ++j)
    points(j) = -1 + r * (1 + cos(M_PI * double(j) / double(m_terms - 1))) / 2;

  for (int j = m_terms; j < 2 * m_terms; ++j)
    points(j) =
        1 - r +
        r * (1 + cos(M_PI * double(j - m_terms) / double(m_terms - 1))) / 2;

  Matrix<double> mat = Matrix<double>::Zero(2 * m_terms, m_terms);
  for (int j = 0; j < m_terms; ++j) {
    Vector<double> values;
    Vector<double> weights = Vector<double>::Zero(m_terms);
    weights(j) = 1;
    evaluateInterpolationPolynomial(weights, points, values);
    mat.col(j) = values;
  }
  return mat;
}
}

#endif
