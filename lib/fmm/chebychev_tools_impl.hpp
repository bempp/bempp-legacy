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

  m_barycentricWeights.resize(m_terms);
  int sign = 1;
  for (int k = 0; k < m_terms; ++k) {
    m_barycentricWeights[k] = sign;
    sign *= -1;
  }
  m_barycentricWeights[0] = double(1) / 2;
  m_barycentricWeights[order] *= double(1) / 2;

  Vector<double> diffWeights = m_barycentricWeights;
  diffWeights[0] *= 4;
  diffWeights[order] *= 4;
  Vector<double> invDiffWeights = Vector<double>::Zero(m_terms);
  invDiffWeights.array() = 1. / diffWeights.array();
  Matrix<double> X = m_nodes.replicate(1, m_terms);
  Matrix<double> dX = X - X.transpose();
  m_chebDiffMatrix = (diffWeights * invDiffWeights.transpose()).eval().array() /
                     (dX + Matrix<double>::Identity(m_terms, m_terms)).array();
  Matrix<double> diag = m_chebDiffMatrix.rowwise().sum().asDiagonal();
  m_chebDiffMatrix.array() -= diag.array();
}

inline const Vector<double> &ChebychevTools::chebychevNodes() const {
  return m_nodes;
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
ChebychevTools::childInterpolationMatrix(double ratio) const {

  double r = 2 * ratio;
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

inline Vector<double>
ChebychevTools::derivativeWeights(const Vector<double> &weights) const {

  return m_chebDiffMatrix * weights;
}

inline Vector<double>
ChebychevTools::derivativeWeights3d(const Vector<double> &weights,
                                    int direction) const {

  Vector<double> result(m_terms * m_terms * m_terms);
  if (direction == 0) {
    for (int i = 0; i < m_terms; ++i)
      for (int j = 0; j < m_terms; ++j) {
        int index = i * m_terms * m_terms + j * m_terms;
        result.segment(index, m_terms) =
            derivativeWeights(weights.segment(index, m_terms));
      }
  } else if (direction == 1) {
    Vector<double> coeffs(m_terms);
    Vector<double> tmp(m_terms);
    for (int i = 0; i < m_terms; ++i) {
      for (int k = 0; k < m_terms; ++k) {
        int index = i * m_terms * m_terms + k;
        for (int j = 0; j < m_terms; ++j)
          coeffs(j) = weights(index + j * m_terms);
        tmp = derivativeWeights(coeffs);
        for (int j = 0; j < m_terms; ++j)
          result(index + j * m_terms) = tmp(j);
      }
    }
  } else if (direction == 2) {

    Vector<double> coeffs(m_terms);
    Vector<double> tmp(m_terms);
    for (int j = 0; j < m_terms; ++j) {
      for (int k = 0; k < m_terms; ++k) {
        int index = j * m_terms + k;
        for (int i = 0; i < m_terms; ++i)
          coeffs(i) = weights(index + i * m_terms * m_terms);
        tmp = derivativeWeights(coeffs);
        for (int i = 0; i < m_terms; ++i)
          result(index + i * m_terms * m_terms) = tmp(i);
      }
    }
  }
  return result;
}
}
#endif
