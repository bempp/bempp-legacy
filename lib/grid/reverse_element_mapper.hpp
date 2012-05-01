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
