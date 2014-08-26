// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_GEOMETRY_INTERFACE_HPP
#define HMAT_GEOMETRY_INTERFACE_HPP

#include "bounding_box.hpp"
#include <array>

namespace hmat {

struct GeometryDataType;

class GeometryInterface {
public:
  virtual shared_ptr<const GeometryDataType> next() const = 0;
  virtual int numberOfEntities() const = 0;
};
}

#endif
