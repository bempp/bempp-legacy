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

#ifndef bempp_mapper_hpp
#define bempp_mapper_hpp

namespace Bempp
{

template <int codim> class Entity;

/** \brief Abstract mapper class.

  A mapper provides a mapping from a set of entities to a range of integers from
  0 to (size() - 1), where size() is the number of entities contained in the set.
  */
class Mapper
{
public:    
    /** \brief Total number of entities in the entity set managed by the mapper. */
    virtual int size() const = 0;

    /** \brief Index of the entity \e of codimension 0.

     The result of calling this method with an entity that does not belong to
     the mapped set is undefined.

     \return An index in the range 0 ... (size() - 1). */
    virtual int entityIndex(const Entity<0>& e) const = 0;
    /** \brief Index of the entity \e of codimension 1.

     \overload
     */
    virtual int entityIndex(const Entity<1>& e) const = 0;
    /** \brief Index of the entity \e of codimension 2.

     \overload
     */
    virtual int entityIndex(const Entity<2>& e) const = 0;
    /** \brief Index of the entity \e of codimension 3.

     \overload
     */
    virtual int entityIndex(const Entity<3>& e) const = 0;
    /** \brief Index of \p i'th subentity of codimension \p codimSub of
     entity \p e of codimension 0. */

    virtual int subEntityIndex(const Entity<0>& e, int i,
                               unsigned int codimSub) const = 0;
};

} // namespace Bempp

#endif
