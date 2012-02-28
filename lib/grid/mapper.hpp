#ifndef MAPPER_HPP
#define MAPPER_HPP

namespace Bempp
{

template <int codim> class Entity;

template <int codim>
class Mapper
{
public:
    /** \brief Map entity to array index. */
    virtual int entityIndex(const Entity<codim>& e) const = 0;
    virtual int subEntityIndex(const Entity<codim>& e, int i,
                               unsigned int codimSub) const = 0;
};

} // namespace Bempp

#endif // MAPPER_HPP
