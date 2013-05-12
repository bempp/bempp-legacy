#ifndef bempp_bounding_box_hpp
#define bempp_bounding_box_hpp

#include "common.hpp"
#include "types.hpp"

namespace Bempp
{

template <typename CoordinateType>
struct BoundingBox
{
    Point3D<CoordinateType> reference;
    Point3D<CoordinateType> lbound;
    Point3D<CoordinateType> ubound;
};

} // namespace Bempp

#endif
