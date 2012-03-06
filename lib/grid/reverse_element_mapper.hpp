#ifndef bempp_reverse_element_mapper_hpp
#define bempp_reverse_element_mapper_hpp

#include <boost/ptr_container/ptr_vector.hpp>
#include "../common/types.hpp"

namespace Bempp
{

template <int codim> class EntityPointer;
class GridView;

/** \brief Mapping from codim-0 entity indices to entity pointers. */
class ReverseElementMapper
{
    template <typename DuneGridView> friend class ConcreteGridView;

private:
    const GridView& m_view;
    boost::ptr_vector<boost::nullable<EntityPointer<0> > > m_cache;

private:
    explicit ReverseElementMapper(const GridView& view);

    void update();

public:
    /** \brief Return an entity pointer referring to the codim-0 entity
      having the given index. */
    const EntityPointer<0>& entityPointer(EntityIndex entityIndex) const;
};

} // namespace Bempp

#endif
