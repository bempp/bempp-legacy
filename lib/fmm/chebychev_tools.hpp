#ifndef bempp_fmm_chebychev_tools_hpp
#define bempp_fmm_chebychev_tools_hpp

#include "fmm_common.hpp"

namespace Fmm {

class ChebychevTools {

public:
  ChebychevTools(int order);

  const Vector<double> &chebychevNodes() const;

  const Matrix<double> &chebychevPolValuesAtNodes() const;

  /** \brief Return the interpolation matrix from a parent to two children boxes
   */
  Matrix<double> interpolateToChildren(double parentLength,
                                       double childLength) const;

  /** \brief Evalute the Chebychevpolynomial using the given weights at the
   * given points. */
  void evaluateInterpolationPolynomial(const Vector<double> &weights,
                                       const Vector<double> &evaluationPoints,
                                       Vector<double> &result) const;

private:
  int m_terms;
  Vector<double> m_nodes;
  Vector<double> m_barycentricWeights;
  Matrix<double> m_chebychevPolValuesAtNodes;
};
}

#include "chebychev_tools_impl.hpp"

#endif
