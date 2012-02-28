#ifndef GEOMETRY_FACTORY_HPP
#define GEOMETRY_FACTORY_HPP

#include <memory>

namespace Bempp
{

class Geometry;

class GeometryFactory
{
public:
    typedef Bempp::Geometry Geometry;

    virtual std::auto_ptr<Geometry> make() const = 0;
};

} // namespace Bempp

#endif // GEOMETRY_FACTORY_HPP
