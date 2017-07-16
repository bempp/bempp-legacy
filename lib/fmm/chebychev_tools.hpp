#ifndef bempp_fmm_chebychev_tools_hpp
#define bempp_fmm_chebychev_tools_hpp

#include "fmm_common.hpp"

namespace Fmm {

class ChebychevTools {

public:
  ChebychevTools(int order);

  const Vector<double> &chebychevNodes() const;

  /** \brief Return the interpolation matrix from a parent to two children boxes
    *        with given ratio of child to parent width. The children can overlap
   * (ratio > .5)
   */
  Matrix<double> childInterpolationMatrix(double ratio) const;

  /** \brief Evalute the Chebychevpolynomial using the given weights at the
   * given points. */
  void evaluateInterpolationPolynomial(const Vector<double> &weights,
                                       const Vector<double> &evaluationPoints,
                                       Vector<double> &result) const;

  /** \brief Compute the derivate of the polynomial with the given weights
       and interpolate on Chebychev points. */
  Vector<double> derivativeWeights(const Vector<double> &weights) const;

  Vector<double> derivativeWeights3d(const Vector<double> &weights,
                                     int direction) const;

private:
  int m_terms;
  Vector<double> m_nodes;
  Vector<double> m_barycentricWeights;
  Matrix<double> m_chebDiffMatrix;
};
}

#include "chebychev_tools_impl.hpp"

#endif
