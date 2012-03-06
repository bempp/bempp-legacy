#include "reverse_element_mapper.hpp"

#include "entity.hpp"
#include "entity_iterator.hpp"
#include "grid_view.hpp"
#include "mapper.hpp"

#include <stdexcept>
#include <memory>

namespace Bempp
{

ReverseElementMapper::ReverseElementMapper(const GridView& view) :
    m_view(view) {
}

void ReverseElementMapper::update() {
    const Mapper& mapper = m_view.elementMapper();
    const int elementCount = m_view.entityCount(0);

    std::auto_ptr<EntityIterator<0> > it = m_view.entityIterator<0>();
    m_cache.clear();
    m_cache.reserve(elementCount);
    for (int i = 0; i < elementCount; ++i)
        m_cache.push_back(0);

    while (!it->finished())
    {
        const Entity<0>& entity = it->entity();
        int index = mapper.entityIndex(entity);
        m_cache.replace(index, it->frozen());
        it->next();
    }
}

const EntityPointer<0>& ReverseElementMapper::entityPointer(
        EntityIndex entityIndex) const {
    return m_cache[entityIndex];
}

} // namespace Bempp
