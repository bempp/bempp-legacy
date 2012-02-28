#ifndef bempp_concrete_geometry_factory_hpp
#define bempp_concrete_geometry_factory_hpp

#include "geometry_factory.hpp"
#include "concrete_geometry.hpp"

namespace Bempp
{

/** \brief Factory able to construct an "empty" geometry wrapping a Dune
  geometry of type DuneGeometry.

  \note For internal use (in integrators from the Fiber module). */

template <typename DuneGeometry>
class ConcreteGeometryFactory : public GeometryFactory
{
    virtual std::auto_ptr<Geometry> make() const {
        return std::auto_ptr<Geometry>(
                    new ConcreteGeometry<DuneGeometry>());
    }
};

} // namespace Bempp

#endif
