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

inline IndexSetType sortIndexSet(const IndexSetType &indexSet,
                                 const Geometry &geometry, int dim) {

  IndexSetType sortedIndexSet(indexSet);

  auto sortFun = [&geometry, dim ](int elem1, int elem2)->bool {

    return (geometry[elem1]->center[dim] < geometry[elem2]->center[dim]);
  };

  std::sort(begin(sortedIndexSet), end(sortedIndexSet), sortFun);

  return sortedIndexSet;
}
}

#endif
