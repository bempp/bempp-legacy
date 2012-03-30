#ifndef bempp_geometry_factory_hpp
#define bempp_geometry_factory_hpp

#include "geometry.hpp"

#include <memory>

namespace Bempp
{

/** \brief Abstract geometry factory. */
class GeometryFactory
{
public:
    typedef Bempp::Geometry Geometry;

    virtual std::auto_ptr<Geometry> make() const = 0;
};

} // namespace Bempp

#endif // GEOMETRY_FACTORY_HPP
