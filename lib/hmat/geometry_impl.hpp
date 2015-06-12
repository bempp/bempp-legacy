// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_GEOMETRY_IMPL_HPP
#define HMAT_GEOMETRY_IMPL_HPP

#include "geometry.hpp"
#include "geometry_data_type.hpp"
#include "geometry_interface.hpp"
#include <algorithm>

namespace hmat {

inline void fillGeometry(Geometry &geometry,
                         GeometryInterface &geometryInterface) {

  geometry.reserve(geometryInterface.numberOfEntities());
  shared_ptr<const GeometryDataType> it;
  while ((it = geometryInterface.next()))
    geometry.push_back(it);
}


}

#endif
