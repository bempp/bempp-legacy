// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_GEOMETRY_DATA_TYPE_IMPL_HPP
#define HMAT_GEOMETRY_DATA_TYPE_IMPL_HPP

#include "geometry_data_type.hpp"

namespace hmat {

inline GeometryDataType::GeometryDataType() : boundingBox(), center() {}

inline GeometryDataType::GeometryDataType(const BoundingBox &boundingBox,
                                          const std::array<double, 3> &center)
    : boundingBox(boundingBox), center(center[0],center[1],center[2]) {}

}

#endif
