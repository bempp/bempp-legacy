#ifndef bempp_bounding_box_hpp
#define bempp_bounding_box_hpp

#include "common.hpp"
#include "types.hpp"

namespace Bempp
{

/** \ingroup common
 *  \brief Bounding box with a reference point. */
template <typename CoordinateType>
struct BoundingBox
{
    /** \brief Reference point. */
    Point3D<CoordinateType> reference;
    /** \brief Lower bound. */
    Point3D<CoordinateType> lbound;
    /** \brief Upper bound. */
    Point3D<CoordinateType> ubound;
};

} // namespace Bempp

#endif
