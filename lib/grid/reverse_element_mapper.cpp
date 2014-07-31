// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

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
    const size_t elementCount = m_view.entityCount(0);

    std::unique_ptr<EntityIterator<0> > it = m_view.entityIterator<0>();
    m_cache.clear();
    m_cache.reserve(elementCount);
    for (size_t i = 0; i < elementCount; ++i)
        m_cache.push_back(0);

    while (!it->finished())
    {
        const Entity<0>& entity = it->entity();
        size_t index = mapper.entityIndex(entity);
        m_cache.replace(index, it->frozen().release());
        it->next();
    }
}

const EntityPointer<0>& ReverseElementMapper::entityPointer(
        EntityIndex entityIndex) const {
    return m_cache[entityIndex];
}

} // namespace Bempp
