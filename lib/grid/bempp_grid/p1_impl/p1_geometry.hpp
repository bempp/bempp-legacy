#ifndef bempp_p1_geometry_hpp
#define bempp_p1_geometry_hpp

#include "../../../common/common.hpp"
#include <dune/geometry/multilineargeometry.hh>
#include <cmath>

namespace BemppGrid {

template <int mydim, int cdim, class>
using P1GridGeometry = Dune::CachedMultiLinearGeometry<double, mydim, cdim>;

}

#endif
