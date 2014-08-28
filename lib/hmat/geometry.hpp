// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_GEOMETRY_HPP
#define HMAT_GEOMETRY_HPP

#include "common.hpp"

#include <array>
#include <vector>

namespace hmat {

class GeometryInterface;
struct GeometryDataType;

typedef std::vector<shared_ptr<const GeometryDataType>> Geometry;

void fillGeometry(Geometry &geometry,
                  GeometryInterface &geometryInterface);

IndexSetType sortIndexSet(const IndexSetType &IndexSet,
                          const Geometry &geometry, int dim);
}

#include "geometry_impl.hpp"

#endif
