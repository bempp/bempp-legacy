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

#include "mass_matrix_container_initialiser.hpp"

#include "../common/lazy.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/types.hpp"
#include "../fiber/scalar_traits.hpp"

#include <armadillo>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType> class BasisData;
template <typename CoordinateType> class Expression;
template <typename CoordinateType> class GeometricalData;

}

namespace Bempp
{

class Grid;
template <size_t codim> class Entity;
template <size_t codim> class EntityPointer;

template <typename ValueType> class MassMatrixContainer;

template <typename BasisFunctionType>
class Space
{
public:
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ComplexType;

    // A grid reference is necessary because e.g. when setting the element
    // variant it is necessary to check whether the element is triangular
    // or quadrilateral. Also, requests for element refinement should probably
    // be made via Space rather than via Grid.
    explicit Space(Grid& grid);

    virtual ~Space();

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

    /** \brief Reference to the basis attached to the specified element. */
    virtual const Fiber::Basis<BasisFunctionType>& basis(const Entity<0>& element) const = 0;

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

    /** \brief Assign global degrees of freedom to local degrees of freedom. */
    virtual void assignDofs() = 0;

    /** \brief True if assignDofs() has been called before, false otherwise. */
    virtual bool dofsAssigned() const = 0;

    /** \brief Number of global degrees of freedom.

        \note Must not be called before asignDofs(). */
    virtual size_t globalDofCount() const = 0;
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

    /** @}
        @name Mass matrix and inverse mass matrix management
        @} */

    template <typename ResultType>
    typename boost::enable_if_c<boost::is_same<ResultType, BasisFunctionType>::value ||
                                boost::is_same<ResultType, ComplexType>::value>::type
    applyMassMatrix(const arma::Col<ResultType>& argument,
                    arma::Col<ResultType>& result) const;
    template <typename ResultType>
    typename boost::enable_if_c<boost::is_same<ResultType, BasisFunctionType>::value ||
                                boost::is_same<ResultType, ComplexType>::value>::type
    applyInverseMassMatrix(const arma::Col<ResultType>& argument,
                           arma::Col<ResultType>& result) const;

    /** @}
        @name Debugging
        @} */
    virtual void dumpClusterIds(const char* fileName,
                                const std::vector<unsigned int>& clusterIds) const = 0;
    /** @} */

protected:
    // void resetMassMatrixContainers() const; // TODO

protected:
    Grid& m_grid;

private:
    void applyMassMatrixBasisFunctionType(
            const arma::Col<BasisFunctionType>& argument,
            arma::Col<BasisFunctionType>& result) const;
    void applyMassMatrixComplexType(
            const arma::Col<ComplexType>& argument,
            arma::Col<ComplexType>& result) const;
    void applyInverseMassMatrixBasisFunctionType(
            const arma::Col<BasisFunctionType>& argument,
            arma::Col<BasisFunctionType>& result) const;
    void applyInverseMassMatrixComplexType(
            const arma::Col<ComplexType>& argument,
            arma::Col<ComplexType>& result) const;

    mutable Lazy<MassMatrixContainer<BasisFunctionType>,
    MassMatrixContainerInitialiser<BasisFunctionType, BasisFunctionType> >
    m_bftMassMatrixContainer;
    mutable Lazy<MassMatrixContainer<ComplexType>,
    MassMatrixContainerInitialiser<BasisFunctionType, ComplexType> >
    m_ctMassMatrixContainer;
};

/** \brief Get pointers to Basis objects corresponding to all elements. */
template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases);

} //namespace Bempp

#endif
