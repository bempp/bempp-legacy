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

#ifndef bempp_space_hpp
#define bempp_space_hpp

#include "../common/types.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/scalar_traits.hpp"
#include <armadillo>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType> class BasisData;
template <typename ValueType> class Expression;
template <typename ValueType> class GeometricalData;

}

namespace Bempp
{

class Grid;
template <int codim> class Entity;
template <int codim> class EntityPointer;

template <typename ValueType>
class Space
{
public:
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;

    // A grid reference is necessary because e.g. when setting the element
    // variant it is necessary to check whether the element is triangular
    // or quadrilateral. Also, requests for element refinement should probably
    // be made via Space rather than via Grid.
    explicit Space(Grid& grid) : m_grid(grid)
    {}

    virtual ~Space()
    {}

    /** @name Attributes
    @{ */

    /** \brief Dimension of the surface on which the functions are defined. */
    virtual int domainDimension() const = 0;
    /** \brief Dimension of the codomain of the functions.

    In other words, number of components of the values of the functions.
    (E.g. H1 space -> 1, H(curl) space on a 2D surface -> 2). */
    virtual int codomainDimension() const = 0;

    /** \brief Reference to the grid on which the functions are defined. */
    const Grid& grid() const { return m_grid; }

    /** \brief Get pointers to Basis objects corresponding to specified elements.

        \param[in]   elements  List of pointers to elements whose bases should
                               be retrieved.
        \param[out]  bases     On exit, a vector whose ith item is a pointer to
                               the Basis object corresponding to the ith item of
                               \p elements. */

    virtual void getBases(const std::vector<const EntityPointer<0>*>& elements,
                          std::vector<const Fiber::Basis<ValueType>*>& bases) const = 0;

    virtual const Fiber::Basis<ValueType>& basis(const Entity<0>& element) const = 0;

    /** \brief Expression returning values of the shape functions of this space. */
    virtual const Fiber::Expression<CoordinateType>& shapeFunctionValueExpression() const = 0;

    /** @}
        @name Element order management
        @{ */
    virtual void setElementVariant(const Entity<0>& element, ElementVariant variant) = 0;
    virtual ElementVariant elementVariant(const Entity<0>& element) const = 0;
    // additional functions for e.g. increasing polynomial order of all elements
    // ...

    /** @}
        @name DOF management
        @{ */

    virtual void assignDofs() = 0;
    virtual bool dofsAssigned() const = 0; // returns a flag that is set to true via assignDofs()

    virtual int globalDofCount() const = 0;
    virtual void globalDofs(const Entity<0>& element,
                            std::vector<GlobalDofIndex>& dofs) const = 0;

    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const = 0;

    // This function is used only by the ACA assembler.
    // For the moment, Point will always be 3D, independently from the
    // actual dimension of the space. Once Ahmed's bemcluster is made dimension-
    // independent, we may come up with a more elegant solution.
    virtual void globalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const = 0;
    /** @} */

protected:
    Grid& m_grid;
};

} //namespace Bempp

#endif
