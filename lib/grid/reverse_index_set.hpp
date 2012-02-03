#ifndef REVERSE_INDEX_SET_HPP
#define REVERSE_INDEX_SET_HPP

#include <boost/ptr_container/ptr_map.hpp>
#include "../common/types.hpp"

namespace Bempp
{

template <int codim> class EntityPointer;
class GridView;

/** \brief Mapping from entity indices to entity pointers. */
class ReverseIndexSet
{
    template <typename DuneGridView> friend class ConcreteGridView;

private:
    typedef boost::ptr_map<EntityIndex, EntityPointer<0> > EntityPointerMap;

    const GridView& m_view;
    EntityPointerMap m_map;

private:
    explicit ReverseIndexSet(const GridView& view);

    void update();

public:
    // if it becomes necessary, other codimensions than 0
    // could also be handled.
    /** \brief Return an entity pointer referring to the entity of codimension 0
      having the given index. */
    const EntityPointer<0>& entityPointerCodim0(
            const EntityIndex& entityIndex) const;
};

} // namespace Bempp

#endif // REVERSE_INDEX_SET_HPP
