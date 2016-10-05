#ifndef bempp_p1_geometry_hpp
#define bempp_p1_geometry_hpp

#include "../../../common/common.hpp"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>
#include <cmath>

namespace BemppGrid {

template <int mydim, int cdim, class>
class P1GridGeometry : public Dune::CachedMultiLinearGeometry<double, mydim, cdim>
{

    typedef Dune::CachedMultiLinearGeometry<double, mydim, cdim> Base;

    public:

        P1GridGeometry(const Dune::GeometryType& type, 
                const std::vector<Dune::FieldVector<double, 3>>& vertices)
            : Base(type, vertices) {}


};

}

#endif
