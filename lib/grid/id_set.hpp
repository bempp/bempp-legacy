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

#ifndef bempp_id_set_hpp
#define bempp_id_set_hpp

#include "../common/common.hpp"

#include <boost/utility/enable_if.hpp>

namespace Bempp
{

/** \cond FORWARD_DECL */
template<int codim> class Entity;
/** \endcond */

/** \ingroup grid
    \brief Abstract wrapper of an id set. */
class IdSet
{
public:
    /** \brief Destructor. */
    virtual ~IdSet() {
    }

    /** \brief Id type.

     \internal Sadly, it is necessary to specify this type uniformly for all grid classes.
     */
    typedef size_t IdType;

    /** \brief Id of the entity \p e of codimension 0. */
    virtual IdType entityId(const Entity<0>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 1. */
    virtual IdType entityId(const Entity<1>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 2. */
    virtual IdType entityId(const Entity<2>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 3. */
    virtual IdType entityId(const Entity<3>& e) const = 0;

    /** \brief Id of \p i'th subentity of codimension \p codimSub of entity \p e of codimension 0.
     */
    virtual IdType subEntityId(const Entity<0>& e, size_t i, int codimSub) const = 0;
};

} // namespace Bempp

#endif // ID_SET_HPP
