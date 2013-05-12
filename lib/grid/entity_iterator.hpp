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

#ifndef bempp_entity_iterator_hpp
#define bempp_entity_iterator_hpp

#include "../common/common.hpp"

#include "entity_pointer.hpp"
#include <memory>

namespace Bempp
{

/** \cond FORWARD_DECL */
template<int codim> class Entity;
/** \endcond */

/** \ingroup grid
    \brief Abstract base class for iterators over entities.

    \param codim Codimension of the entities iterated over.

    Typical usage (here \p ConcreteEntityIterator is some subclass of <tt>EntityIterator<codim></tt>):

    \code
    ConcreteEntityIterator* it = ...();
    while (!it->finished())
    {
    const ConcreteEntity& e = it.entity();
    do_stuff(e);
    it->next();
    }
    delete it;
    \endcode

    \internal Reminder: The template parameter codim is necessary because the
    entity() method must know the codimension of the entity to which it returns
    a reference.
*/
template<int codim>
class EntityIterator: public EntityPointer<codim>
{
protected: /* Can't be changed to private, derived classes use it */
    bool m_finished;
public:
    /** \brief Increment iterator. */
    virtual void next() = 0;

    /** \brief True if iterator points past the end of the iteration range */
    bool finished() const {
        return m_finished;
    }

    // virtual void reset() = 0; // Would such a method be useful?

    // This redeclaration appears so that the docstring can be changed wrt. to the base class.
    /** \brief Read-only access to the entity referenced by the iterator. */
    virtual const Entity<codim>& entity() const = 0;

    /** \brief A stable pointer to the entity currently referenced by the iterator.

    The returned pointer is guaranteed to keep referring to the same entity
    even if the iterator is incremented or destroyed, as long as the grid is
    not adapted and the grid object itself stays alive. */
    virtual std::auto_ptr<EntityPointer<codim> > frozen() const = 0;
};

} // namespace Bempp

#endif
