#ifndef bempp_fmm_chebychev_tools_hpp
#define bempp_fmm_chebychev_tools_hpp

#include "fmm_common.hpp"

namespace Fmm {

class ChebychevTools {

public:
  ChebychevTools(int order);

  void chebychevNodes(Vector<double> &nodes) const;

private:
  int m_order;
};
}

#include "chebychev_tools_impl.hpp"

#endif
