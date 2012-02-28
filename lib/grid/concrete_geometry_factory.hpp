#ifndef CONCRETE_GEOMETRY_FACTORY_HPP
#define CONCRETE_GEOMETRY_FACTORY_HPP

#include "geometry_factory.hpp"
#include "concrete_geometry.hpp"

namespace Bempp
{

template <typename DuneGeometry>
class ConcreteGeometryFactory : public GeometryFactory
{
    virtual std::auto_ptr<Geometry> make() const {
        return std::auto_ptr<Geometry>(
                    new ConcreteGeometry<DuneGeometry>());
    }
};

}

#endif // CONCRETE_GEOMETRY_FACTORY_HPP
