// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_GEOMETRY_DATA_TYPE_HPP
#define HMAT_GEOMETRY_DATA_TYPE_HPP

#include "common.hpp"
#include "geometry_interface.hpp"
#include "bounding_box.hpp"

#include <array>
#include <vector>

namespace hmat {
struct GeometryDataType {

  GeometryDataType();
  GeometryDataType(const BoundingBox &boundingBox,
                   const std::array<double, 3> &center);

  BoundingBox boundingBox;
  Point center;
};
}

#include "geometry_data_type_impl.hpp"

#endif
