#ifndef bempp_fmm_chebychev_tools_impl_hpp
#define bempp_fmm_chebychev_tools_impl_hpp

#include "chebychev_tools.hpp"

#include <cmath>

namespace Fmm {

inline ChebychevTools::ChebychevTools(int order) : m_order(order) {}

inline void ChebychevTools::chebychevNodes(Vector<double> &nodes) const {

  nodes.resize(m_order);
  for (int m = 0; m < m_order; ++m)
    nodes(m) = cos(M_PI * double(2 * m + 1) / double(2 * m_order));
}
}

#endif
