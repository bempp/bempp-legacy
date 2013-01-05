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

#include "../common/common.hpp"
#include "bempp/common/config_trilinos.hpp"

#include "../common/deprecated.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"
#include "../fiber/scalar_traits.hpp"

#include "../common/armadillo_fwd.hpp"
#include <vector>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ValueType> class Basis;
template <typename ValueType> class BasisData;
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename CoordinateType> class GeometricalData;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class Grid;
template <int codim> class Entity;
template <int codim> class EntityPointer;
template <typename ValueType> class DiscreteBoundaryOperator;
/** \endcond */

enum DofType
{
    GLOBAL_DOFS, FLAT_LOCAL_DOFS
};

/** \ingroup space
 *  \brief Function space.
 *
 *  This class represents a space of functions defined on a grid. The space is
 *  spanned by a finite number of scalar- or vector-valued basis functions. The
 *  template parameter \p BasisFunctionType is the type of the values of (the
 *  components of) these basis functions, and can be set to \c float,
 *  \c double, <tt>std::complex<float></tt> or <tt>std::complex<double></tt>.
 *
 *  The basis functions of a space, also known as global degrees of freedom
 *  (DOFs), can have support extending over multiple elements of the grid. They
 *  are, however, composed of one or more *local* basis functions (local
 *  degrees of freedom), each of which resides on a single element. The mapping
 *  of local to global degrees of freedom is triggered by calling the function
 *  assignDofs(). Many other member functions of Space may only be invoked after
 *  assignDofs() has beeen called. */
template <typename BasisFunctionType>
class Space
{
public:
    /** \brief Type used to represent coordinates. */
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    /** \brief Equivalent to std::complex<CoordinateType>. */
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ComplexType;
    /** \brief Appropriate instantiation of Fiber::CollectionOfBasisTransformations. */
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;

    // A grid reference is necessary because e.g. when setting the element
    // variant it is necessary to check whether the element is triangular
    // or quadrilateral. Also, requests for element refinement should probably
    // be made via Space rather than via Grid.
    /** \brief Constructor.
     *
     *  \param[in] grid Grid on which functions from this space should be
     *  defined. */
    explicit Space(const shared_ptr<const Grid>& grid);

    /** \brief Destructor. */
    virtual ~Space();

    /** @name Attributes
    @{ */

    /** \brief Dimension of the grid on which functions from this space are
     *  defined. */
    virtual int domainDimension() const = 0;

    /** \brief Dimension of the codomain of the functions.
     *
     * In other words, number of components of the values of the functions.
     * (E.g. H1 space -> 1, H(curl) space on a 2D surface -> 2). */
    virtual int codomainDimension() const = 0;

    /** \brief Shared pointer to the grid on which the functions from this space
     *  are defined. */
    shared_ptr<const Grid> grid() const { return m_grid; }

    /** \brief Reference to the basis attached to the specified element. */
    virtual const Fiber::Basis<BasisFunctionType>& basis(
            const Entity<0>& element) const = 0;

    /** \brief Transformation mapping basis functions to shape functions.
     *
     *  This function returns a CollectionOfBasisTransformations object
     *  consisting of a single transformation that maps values of basis
     *  functions defined on a reference element to those of *shape functions*
     *  defined on a particular element of the grid.
     *
     *  This transformation is the identity for spaces of scalar-valued
     *  functions, but may be more complicated for spaces of vector-valued
     *  functions, e.g. \f$H(\mathrm{curl})\f$.
     *
     *  \todo Perhaps change the name of this method to something more
     *  understandable, like basisToShapeFunctionTransformation. */
    virtual const CollectionOfBasisTransformations&
    shapeFunctionValue() const = 0;

    /** @}
        @name Element order management
        @{ */

    /** \brief Set the variant of element \p element to \p variant.
     *
     *  The element variant determines the set of basis functions defined on
     *  the element (e.g. maximum polynomial order). Different subclasses of
     *  Space interpret the \p variant argument in different ways; for more
     *  information, see the documentation of these subclasses.
     *
     *  \note Calling this function only makes sense for subclasses implementing
     *  adaptive function spaces. Currently there are no such subclasses. */
    virtual void setElementVariant(const Entity<0>& element, ElementVariant variant) = 0;

    /** \brief Return current variant of element \p element.
     *
     *  See the documentation of setElementVariant() for more information. */
    virtual ElementVariant elementVariant(const Entity<0>& element) const = 0;

    // additional functions for e.g. increasing polynomial order of all elements
    // ...

    /** @}
        @name DOF management
        @{ */

    /** \brief Assign global degrees of freedom to local degrees of freedom.
     *
     *  \deprecated It is not necessary any more to call this function, since
     *  DOF assignment is now done automatically in the constructor.
     */
    BEMPP_DEPRECATED void assignDofs();

    /** \brief True if global degrees of freedom have been already assigned to
     *  local degrees of freedom, false otherwise.
     *
     *  \deprecated Since DOF assignment is now done automatically in the
     *  constructor, this function always returns true.
     */
    BEMPP_DEPRECATED bool dofsAssigned() const;

    /** \brief Total number of local degrees of freedom on all elements. */
    virtual size_t flatLocalDofCount() const = 0;

    /** \brief Number of global degrees of freedom. */
    virtual size_t globalDofCount() const = 0;

    /** \brief Map local degrees of freedom residing on an element to global
     *  degrees of freedom.
     *
     *  \param[in] element An element of the grid grid().
     *  \param[out] dofs   Indices of the global degrees of freedom
     *                     corresponding to the local degrees of freedom
     *                     residing on \p element. */
    virtual void getGlobalDofs(const Entity<0>& element,
                               std::vector<GlobalDofIndex>& dofs) const;

    virtual void getGlobalDofs(const Entity<0>& element,
                               std::vector<GlobalDofIndex>& dofs,
                               std::vector<BasisFunctionType>& localDofWeights) const;

    /** \brief Map global degrees of freedom to local degrees of freedom.
     *
     *  \param[in] globalDofs
     *     Vector containing indices of global degrees of freedom.
     *
     *  \param[out] localDofs
     *     Vector whose <tt>i</tt>th element is the vector containing all the
     *     local degrees of freedom that are mapped to the global
     *     degree of freedom <tt>globalDofs[i]</tt>.
     *
     *  Note that a local degree of freedom (LocalDof) is a combination of an
     *  EntityIndex and LocalDofIndex, as explained in its documentation. */
    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const;

    virtual void global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs,
            std::vector<std::vector<BasisFunctionType> >& localDofWeights) const;

    /** \brief Map flat indices of local degrees of freedom to local degrees of freedom.
     *
     *  \param[in] flatLocalDofs
     *     Vector containing flat indices of local degrees of freedom.
     *
     *  \param[out] localDofs
     *     Vector whose <tt>i</tt>th element is the local degree of freedom
     *     with flat index given by <tt>flatLocalDofs[i]</tt>. */
    virtual void flatLocal2localDofs(
            const std::vector<FlatLocalDofIndex>& flatLocalDofs,
            std::vector<LocalDof>& localDofs) const = 0;

    // These functions are used only by the ACA assembler.
    // For the moment, Point will always be 3D, independently from the
    // actual dimension of the space. Once Ahmed's bemcluster is made dimension-
    // independent, we may come up with a more elegant solution.
    /** \brief Retrieve positions of global degrees of freedom.
     *
     *  \param[out] positions
     *    Vector whose <em>i</em>th element contains the coordinates
     *    of the point taken to be the "position" (in some sense) of
     *    <em>i</em>th global degree of freedom.
     *
     *  \note This function is intended as a helper for clustering algorithms
     *  used in matrix compression algorithms such as adaptive cross
     *  approximation. */
    virtual void getGlobalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const = 0;

    /** \brief Retrieve positions of local degrees of freedom ordered by their
     *  flat index.
     *
     *  \param[out] positions
     *    Vector whose <em>i</em>th element contains the coordinates
     *    of the point taken to be the ``position'' (in some sense) of
     *    the local degree of freedom with flat index <em>i</em>.
     *
     *  \note This function is intended as a helper for clustering algorithms
     *  used in matrix compression algorithms such as adaptive cross
     *  approximation. */
    virtual void getFlatLocalDofPositions(
            std::vector<Point3D<CoordinateType> >& positions) const = 0;

    /** @}
        @name Debugging
        @} */

    /** \brief Write a VTK file showing the distribution of global or
     *  flat local degrees of freedom into clusters.
     *
     *  \param[in] fileName
     *    Name of the VTK file to be created (without extension).
     *  \param[in] clusterIdsOfGlobalDofs
     *    Vector whose <em>i</em>th element contains the identifier of the
     *    cluster to which <em>i</em>th global degree of freedom has been assigned.
     *
     *  This function generates a VTK file containing a single data series
     *  mapping the ``positions'' (see getGlobalDofPositions()) of global degrees
     *  of freedom to the identifiers of the clusters to which these degrees of
     *  freedom have been assigned. It is intended for debugging clustering
     *  algorithms.
     *
     *  \deprecated This function is deprecated. Use dumpClusterIdEx()
     *  instead, which supports dumping of cluster identifiers of flat
     *  local degrees of freedom in addition to the global ones.
     */
    BEMPP_DEPRECATED virtual void dumpClusterIds(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs) const = 0;

    /** \brief Write a VTK file showing the distribution of global or
     *  flat local degrees of freedom into clusters.
     *
     *  \param[in] fileName
     *    Name of the VTK file to be created (without extension).
     *  \param[in] clusterIdsOfGlobalDofs
     *    Vector whose <em>i</em>th element contains the identifier of the
     *    cluster to which <em>i</em>th degree of freedom has been assigned.
     *  \param[in] dofType
     *    Type of degrees of freedom (GLOBAL_DOFS or FLAT_LOCAL_DOFS).
     *
     *  This function generates a VTK file containing a single data
     *  series mapping the ``positions'' (see getGlobalDofPositions()
     *  and getFlatLocalDofPositions()) of the chosen type of degrees
     *  of freedom to the identifiers of the clusters to which these
     *  degrees of freedom have been assigned. It is intended for
     *  debugging clustering algorithms.
     *
     *  This function supersedes dumpClusterIds().
     */
    virtual void dumpClusterIdsEx(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
            DofType dofType) const;
    /** @} */

private:
    /** \cond PRIVATE */
    shared_ptr<const Grid> m_grid;
    /** \endcond */
};

/** \brief Get pointers to Basis objects corresponding to all elements of the grid
 *  on which a function space is defined.
 *
 *  \param[in] space
 *    A Space object.
 *  \param[out] bases
 *    Vector whose <em>i</em>th element is a pointer to the Basis object
 *    representing the local degrees of freedom residing on the <em>i</em>th
 *    element of the grid on which \p space is defined.
 *
 *  An exception is raised if <tt>space.assignDofs()</tt> has not been called
 *  prior to calling this function. */
template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases);

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
constructOperatorMappingGlobalToFlatLocalDofs(const Space<BasisFunctionType>& space);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
constructOperatorMappingFlatLocalToGlobalDofs(const Space<BasisFunctionType>& space);
#endif // WITH_TRILINOS

} //namespace Bempp

#endif
