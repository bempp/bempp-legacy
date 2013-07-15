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

#ifndef bempp_entity_hpp
#define bempp_entity_hpp

#include "../common/common.hpp"

#include "geometry_type.hpp"

#include <boost/utility/enable_if.hpp>
#include <dune/common/static_assert.hh>
#include <dune/grid/common/geometry.hh>
#include <memory>
#include <stdexcept>

namespace Bempp
{

/** \cond FORWARD_DECL */
class Geometry;
template<int codim> class EntityPointer;
template<int codim> class EntityIterator;
/** \endcond */

/** \ingroup grid
    \brief Abstract wrapper of an entity of codimension \p codim. */
template<int codim>
class Entity
{
public:
    /** \brief Destructor */
    virtual ~Entity() {
    }

    /** \brief %Entity codimension. */
    enum {
        codimension = codim
    };

    /** \brief %Entity level. */
    virtual size_t level() const = 0;

    /** \brief Reference to the geometry of this entity.

     This object gives, among other things, the map from a reference element to world coordinates.

     \note Be careful when storing such references. If the state
     of any object is changed, e.g. an iterator is advanced, there
     is no guarantee that the reference remains valid.
     */
    virtual const Geometry& geometry() const = 0;

    /** \brief Type of the reference element.

     The type can be used to access the Dune::GenericReferenceElement.
     */
    virtual GeometryType type() const = 0;
};

/** \ingroup grid
    \brief Abstract wrapper of an entity of codimension 0. */
template<>
class Entity<0>
{
public:
    /** \brief Destructor. */
    virtual ~Entity() {
    }

    /** \brief %Entity codimension. */
    enum {
        codimension = 0
    };

    /** @name Methods shared by entities of all codimensions
    @{ */

    /** \brief %Entity level */
    virtual size_t level() const = 0;

    /** \brief Reference to the geometry of this entity.

     This object gives, among other things, the map from a reference element to world coordinates.

     \note Be careful when storing such references. If the state
     of any object is changed, e.g. an iterator is advanced, there
     is no guarantee that the reference remains valid.
     */
    virtual const Geometry& geometry() const = 0;

    /** \brief Type of the reference element.

    The type can be used to access the Dune::GenericReferenceElement.
     */
    virtual GeometryType type() const = 0;

    /** @}
    @name Extended interface of entities of codimension 0
    @{ */

    /** \brief Number of subentities of codimension \p codimSub. */
    // Default implementation, specialisations for potentially allowed
    // codimensions (1 to 3) follow after class declaration.
    template<int codimSub> size_t subEntityCount() const {
        return 0;
    }

    /** \brief Iterator over subentities of codimension \p codimSub.

    \note codimSub must be greater than 0 and less than the dimension of the grid. */
    // Default implementation, specialisations for potentially allowed codimensions follow
    // after class declaration.
    template<int codimSub>
    std::auto_ptr<EntityIterator<codimSub> > subEntityIterator() const {
        throw std::logic_error("Entity::subEntityIterator(): invalid entity codimension");
    }

    /** \brief Inter-level access to father entity on the next-coarser grid.

     The given entity resulted directly from a subdivision of its father
     entity. For macro (level-0) elements a null auto_ptr is returned.
     */
    virtual std::auto_ptr<EntityPointer<0> > father() const = 0;

    /** \brief True if entity has a father entity which can be accessed
     using the father() method.
     */
    virtual bool hasFather() const = 0;

    /** \brief True if the entity is contained in the leaf grid. */
    virtual bool isLeaf() const = 0;

    /** \brief True if the element is of regular type in red/green type refinement.

     In bisection or hanging node refinement this is always true.
     */
    virtual bool isRegular() const = 0;

    /** \brief Inter-level access to elements that resulted from (recursive)
     subdivision of this element.

     \param[in] maxlevel Iterator does not stop at elements with level greater than \p maxlevel.
     \return Iterator to the first son (level is not greater than \p maxlevel)
     */
    virtual std::auto_ptr<EntityIterator<0> > sonIterator(int maxlevel) const = 0;

    /** \brief True if the entity has been created during the last call to adapt().
     */
    virtual bool isNew() const = 0;

    /** \brief True if the entity might disappear during the next call to adapt().
     *
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    virtual bool mightVanish() const = 0;

    // Deferred for possible later implementation:
    // * Iteration over neighbours: Dune methods ileafbegin(), ileafend(), ilevelbegin(), ilevelend()
    // * Information about the way this element has been subdivided from
    //   its father element: Dune method geometryInFather().

private:
    /** @}
    @name Internal methods
    @{ */

    /** \brief Iterator over subentities of codimension 1. */
    virtual std::auto_ptr<EntityIterator<1> > subEntityCodim1Iterator() const = 0;
    /** \brief Iterator over subentities of codimension 2. */
    virtual std::auto_ptr<EntityIterator<2> > subEntityCodim2Iterator() const = 0;
    /** \brief Iterator over subentities of codimension 3. */
    virtual std::auto_ptr<EntityIterator<3> > subEntityCodim3Iterator() const = 0;

    /** \brief Number of subentities of codimension 1. */
    virtual size_t subEntityCodim1Count() const = 0;
    /** \brief Number of subentities of codimension 1. */
    virtual size_t subEntityCodim2Count() const = 0;
    /** \brief Number of subentities of codimension 1. */
    virtual size_t subEntityCodim3Count() const = 0;

    /** @} */
};

} // namespace Bempp

#include "entity_iterator.hpp"

namespace Bempp
{

template<>
inline std::auto_ptr<EntityIterator<1> > Entity<0>::subEntityIterator<1>() const
{
    return subEntityCodim1Iterator();
}
template<>
inline std::auto_ptr<EntityIterator<2> > Entity<0>::subEntityIterator<2>() const
{
    return subEntityCodim2Iterator();
}
template<>
inline std::auto_ptr<EntityIterator<3> > Entity<0>::subEntityIterator<3>() const
{
    return subEntityCodim3Iterator();
}

template<>
inline size_t Entity<0>::subEntityCount<1>() const
{
    return subEntityCodim1Count();
}
template<>
inline size_t Entity<0>::subEntityCount<2>() const
{
    return subEntityCodim2Count();
}
template<>
inline size_t Entity<0>::subEntityCount<3>() const
{
    return subEntityCodim3Count();
}

} // namespace Bempp

#endif
