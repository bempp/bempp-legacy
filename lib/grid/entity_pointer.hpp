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
#ifndef bempp_entity_pointer_hpp
#define bempp_entity_pointer_hpp

#include "../common/common.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
template<int codim> class Entity;
/** \endcond */

/**
 \brief Abstract base class for an object providing read-only access to an
 entity of codimension \p codim.

 \internal Reminder: The template parameter codim is necessary because the
 entity() method must know the codimension of the entity to which it returns
 a reference.
 */
template<int codim>
class EntityPointer
{
public:
    /** \brief Destructor */
    virtual ~EntityPointer() {
    }

    /** \brief Entity codimension */
    enum {
        codimension = codim
    };

    /** \brief Read-only access to the underlying entity */
    virtual const Entity<codim>& entity() const = 0;
};

} // namespace Bempp

#endif
