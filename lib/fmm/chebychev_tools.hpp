#ifndef bempp_fmm_chebychev_tools_hpp
#define bempp_fmm_chebychev_tools_hpp

#include "fmm_common.hpp"

namespace Fmm {

class ChebychevTools {

public:
  ChebychevTools(int order);

  const Vector<double> &chebychevNodes() const;

  const Matrix<double> &chebychevPolValuesAtNodes() const;

private:
  int m_terms;
  Vector<double> m_nodes;
  Matrix<double> m_chebychevPolValuesAtNodes;
};
}

#include "chebychev_tools_impl.hpp"

#endif
