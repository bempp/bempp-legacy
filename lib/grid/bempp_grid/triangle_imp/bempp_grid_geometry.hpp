#ifndef bempp_grid_triangle_impl_geometry_hpp
#define bempp_grid_triangle_impl_geometry_hpp

#include "../../../common/common.hpp"
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>
#include <cmath>

namespace BemppGrid {

template <int mydim, int cdim, class>
class Geometry : public Dune::MultiLinearGeometry<double, mydim, cdim>
{

    typedef Dune::MultiLinearGeometry<double, mydim, cdim> Base;

    public:

        Geometry(const Dune::GeometryType& type, 
                const std::vector<Dune::FieldVector<double, 3>>& vertices)
            : Base(type, vertices) {}


};

}

#endif
