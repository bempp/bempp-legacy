#include "reverse_index_set.hpp"

#include "entity.hpp"
#include "entity_iterator.hpp"
#include "grid_view.hpp"
#include "index_set.hpp"

#include <stdexcept>
#include <memory>

namespace Bempp
{

ReverseIndexSet::ReverseIndexSet(const GridView& view) :
    m_view(view) {
}

void ReverseIndexSet::update() {
    const IndexSet& indexSet = m_view.indexSet();
    std::auto_ptr<EntityIterator<0> > it = m_view.entityIterator<0>();

    GeometryType type;
    int index;
    while (!it->finished())
    {
        const Entity<0>& entity = it->entity();
        type = entity.type();
        index = indexSet.entityIndex(entity);
        m_map.insert(EntityIndex(type, index), it->frozen());
        it->next();
    }
}

const EntityPointer<0>& ReverseIndexSet::entityPointerCodim0(
        const EntityIndex& entityIndex) const {
    return m_map.at(entityIndex);
}

} // namespace Bempp
