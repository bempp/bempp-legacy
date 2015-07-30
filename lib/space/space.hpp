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

#include "space_identifier.hpp"

#include "../common/common.hpp"

#include "../common/bounding_box.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/deprecated.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/collection_of_basis_transformations.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/global_parameters.hpp"
#include "../hmat/cluster_tree.hpp"


#include <vector>


/** \cond FORWARD_DECL */
namespace Fiber {
template <typename ValueType> class BasisData;
template <typename CoordinateType> class GeometricalData;
} // namespace Fiber


/** \endcond */


namespace Bempp {

/** \cond FORWARD_DECL */
class Grid;
class GridView;
class GeometryFactory;
template <int codim> class Entity;
template <int codim> class EntityPointer;
template <typename ValueType> class DiscreteSparseBoundaryOperator;
template <typename ValueType> class DiscreteBoundaryOperator;
/** \endcond */

enum DofType { GLOBAL_DOFS, FLAT_LOCAL_DOFS };

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
template <typename BASISFUNCTIONTYPE> class Space {
public:
  /** \brief Instantiation type **/
  typedef BASISFUNCTIONTYPE BasisFunctionType;
  /** \brief Type used to represent coordinates. */
  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
  /** \brief Equivalent to std::complex<CoordinateType>. */
  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ComplexType;
  /** \brief Appropriate instantiation of
   * Fiber::CollectionOfShapesetTransformations. */
  typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
      CollectionOfShapesetTransformations;
  /** \brief Appropriate instantiation of
   * Fiber::CollectionOfBasisTransformations. */
  typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
      CollectionOfBasisTransformations;

  // A grid reference is necessary because e.g. when setting the element
  // variant it is necessary to check whether the element is triangular
  // or quadrilateral. Also, requests for element refinement should probably
  // be made via Space rather than via Grid.
  /** \brief Constructor.
   *
   *  \param[in] grid Grid on which functions from this space should be
   *  defined.
   *
   *  \param[in] level Level of the grid on which the space should be defined.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  explicit Space(const shared_ptr<const Grid> &grid);

  /** \brief Copy Constructor */
  Space(const Space<BasisFunctionType> &other);

  /** \brief Destructor. */
  virtual ~Space();

  /** \brief Assignment operator */
  Space<BasisFunctionType> &operator=(const Space<BasisFunctionType> &other);

  /** @name Attributes
  @{ */

  /** \brief Return a shared pointer to an appropriate counterpart to this
   *  space, with basis functions extending only over single elements.
   *
   *  Let \f$S = span_{n=1}^N f_n\f$ denote the function space represented by
   *  this object, with \f$f_n\f$ its basis functions.
   *
   *  If the functions \f$f_n\f$ are scalar-valued, discontinuousSpace()
   *  should return a shared pointer to an object representing a space \f$T =
   *  span_{m=1}^M g_m\f$ with basis functions \f$g_m\f$ such that:
   *
   *  1. \f$T = S\f$.
   *
   *  2. The support of each basis function \f$g_m\f$ is a single element.
   *
   *  3. Each function \f$f_n\f$ has a unique representation in the basis of
   *     \f$\{g_m\}_{m=1}^M\f$ and each function \f$g_m\f$ contributes to
   *exactly
   *     one function \f$f_n\f$.
   *
   *  If the values of functions \f$f_n\f$ are vectors with \f$d\f$
   *  components, discontinuousSpace() should return a shared pointer to an
   *  object representing a space \f$T = span_{m=1}^M g_m\f$ of with
   *  *scalar-valued* basis functions \f$g_m\f$ such that
   *
   *  1. For each \f$i = 1, 2, \cdots, d\f$ it holds that \f$T \supset
   *     span_{n=1}^N (f_n)_i\f$, where \f$(f_n)_i\f$ denotes the \f$i\f$th
   *     component of \f$f_n\f$.
   *
   *  2. The support of each basis function \f$g_m\f$ is a single element.
   *
   *  \param[in] self This must be a shared pointer to <tt>*this</tt>.
   */
  virtual shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const = 0;

  /** \brief Return true if each basis function of this space extends over
   *  only a single element, false otherwise. */
  virtual bool isDiscontinuous() const = 0;

  /** \brief Return \p true if space is based on a barycentric refinement */

  virtual bool isBarycentric() const = 0;

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
  virtual shared_ptr<const Grid> grid() const { return m_grid; }

  /** \brief Reference to the shapeset attached to the specified element.
   *
   *  \deprecated This function is deprecated and will be removed in a future
   *  release of BEM++. Use shapeset() instead.
   */
  virtual BEMPP_DEPRECATED const Fiber::Basis<BasisFunctionType> &
  basis(const Entity<0> &element) const {
    // It might be good to print a deprecation warning
    return dynamic_cast<const Fiber::Basis<BasisFunctionType> &>(
        shapeset(element));
  }

  /** \brief Reference to the shapeset attached to the specified element.
   */
  virtual const Fiber::Shapeset<BasisFunctionType> &
  shapeset(const Entity<0> &element) const {
    throw NotImplementedError(
        "Space::shapeset(): not implemented.\nNote that the "
        "Space::basis() function has been renamed to shapeset(). "
        "If you have implemented basis() in a subclass of Space, "
        "please implement shapeset() instead.");
  }

  /** \brief Transformation mapping shape functions to basis functions.
   *
   *  \deprecated This function is deprecated and will be removed in a future
   *  release of BEM++. Use basisFunctionValue() instead.
   */
  virtual BEMPP_DEPRECATED const CollectionOfBasisTransformations &
  shapeFunctionValue() const {
    return dynamic_cast<const CollectionOfBasisTransformations &>(
        basisFunctionValue());
  }

  /** \brief Return an equivalent space (in terms of global Dofs), but defined
   * using
   *  local dofs on the barycentrically refined grid. */
  virtual shared_ptr<const Space<BasisFunctionType>> barycentricSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;

  /** \brief Return the grid level of the current space */
  virtual unsigned int level() const { return m_level; }

  /** \brief Return the underlying grid dimension */
  virtual int gridDimension() const;

  /** \brief Return the underlying world dimension */
  virtual int worldDimension() const;

  /** \brief Return the grid view of the current space */
  virtual const GridView &gridView() const;

  /** \brief Transformation mapping shape functions to basis functions.
   *
   *  This function returns a CollectionOfShapesetTransformations object
   *  consisting of a single transformation that maps values of shape
   *  functions defined on a reference element to those of *basis functions*
   *  defined on a particular element of the grid.
   *
   *  This transformation is the identity for spaces of scalar-valued
   *  functions, but may be more complicated for spaces of vector-valued
   *  functions, e.g. \f$H(\mathrm{curl})\f$.
   *
   */
  virtual const CollectionOfShapesetTransformations &
  basisFunctionValue() const {
    throw NotImplementedError(
        "Space::basisFunctionValue(): not implemented.\n"
        "Note that the Space::shapeFunctionValue() function has "
        "been renamed to basisFunctionValue(). If you have "
        "implemented shapeFunctionValue() in a subclass of Space, "
        "please implement basisFunctionValue() instead.");
  }

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
  virtual void setElementVariant(const Entity<0> &element,
                                 ElementVariant variant) = 0;

  /** \brief Return current variant of element \p element.
   *
   *  See the documentation of setElementVariant() for more information. */
  virtual ElementVariant elementVariant(const Entity<0> &element) const = 0;

  /** \brief Return the GeometryFactory associated with the mesh. */

  virtual shared_ptr<GeometryFactory> elementGeometryFactory() const {
    return m_elementGeometryFactory;
  }

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
   *  \param[in] element
   *    An element of the grid grid().
   *  \param[out] dofs
   *    Vector whose <em>i</em>th element is the index of the global degrees
   *    of freedom to which the <em>i</em>th local degree of freedom residing
   *    on \p element contributes. A negative number means that a given local
   *    degree of freedom does not contribute to any global one.
   *
   *  \note This function is deprecated. Use the other overload taking the
   *  additional output parameter \p localDofWeights.
   */
  virtual void getGlobalDofs(const Entity<0> &element,
                             std::vector<GlobalDofIndex> &dofs) const;

  /** \brief Map local degrees of freedom residing on an element to global
   *  degrees of freedom.
   *
   *  \param[in] element
   *    An element of the grid grid().
   *  \param[out] dofs
   *    Vector whose <em>i</em>th element is the index of the global degrees
   *    of freedom to which the <em>i</em>th local degree of freedom residing
   *    on \p element contributes. A negative number means that a given local
   *    degree of freedom does not contribute to any global one.
   *  \param[out] localDofWeights
   *    Vector whose <em>i</em>th element is the weight with which the
   *    <em>i</em>th local degree of freedom residing on \p element
   *    contributes to "its" global degree of freedom. */
  virtual void
  getGlobalDofs(const Entity<0> &element, std::vector<GlobalDofIndex> &dofs,
                std::vector<BasisFunctionType> &localDofWeights) const;

  /** \brief Return true if both spaces act on the same grid. */
  virtual bool gridIsIdentical(const Space<BasisFunctionType> &other) const;

  /** \brief Return the identifier of the space. */
  virtual SpaceIdentifier spaceIdentifier() const = 0;

  /** \brief Return true if \p other is compatible to this space, i.e. the
   * global
   * dofs of the two spaces agree with each other. */
  virtual bool
  spaceIsCompatible(const Space<BasisFunctionType> &other) const = 0;

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
  virtual void
  global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                   std::vector<std::vector<LocalDof>> &localDofs) const;

  virtual void global2localDofs(
      const std::vector<GlobalDofIndex> &globalDofs,
      std::vector<std::vector<LocalDof>> &localDofs,
      std::vector<std::vector<BasisFunctionType>> &localDofWeights) const;

  /** \brief Map flat indices of local degrees of freedom to local degrees of
   *freedom.
   *
   *  \param[in] flatLocalDofs
   *     Vector containing flat indices of local degrees of freedom.
   *
   *  \param[out] localDofs
   *     Vector whose <tt>i</tt>th element is the local degree of freedom
   *     with flat index given by <tt>flatLocalDofs[i]</tt>. */
  virtual void
  flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                      std::vector<LocalDof> &localDofs) const = 0;

  /** @}
      @name Function interpolation
      @} */

  /** \brief Retrieve the interpolation points of the global degrees of freedom.
   *
   *  This function, along with getNormalsAtGlobalDofInterpolationPoints()
   *  and getGlobalDofInterpolationDirections(), can be used in the
   *  interpolation of functions in finite-dimensional function spaces
   *  represented with classes derived from Space. Let \f$V\f$ be a function
   *  space defined on grid \f$\Gamma\f$ with a basis \f$(f_i)_{i=1}^N\f$ of
   *  functions \f$f_i : \Gamma \to \mathbb{R}^n\f$ (or \f$\mathbb{C}^n\f$.
   *  We say that this basis is *interpolatory* if there exist points \f$x_j
   *  \in \Gamma\f$ and \f$n\f$-dimensional vectors \f$d_j\f$ (\f$j = 1, 2,
   *  \dots, N\f$ such that
   *
   *  \f[
   *      f_i(x_j) \cdot d_j =
   *      \begin{cases}
   *         1 &\text{for} i = j,\\
   *         0 &\text{otherwise}.
   *      \end{cases}
   *  \f]
   *
   *  For any appropriate function \f$u\f$ defined on \f$Gamma\f$ the numbers
   *  \f$u(x_i) \cdot d_i\f$ can then be taken as the expansion coefficients
   *  in the basis \f$(f_i)_{i=1}^N\f$ of its interpolation. ("Appropriate"
   *  means that if, for example, the space \f$V\f$ is a space of functions
   *  tangential to the grid, then \f$u\f$ should also be tangential to the
   *  grid.)
   *
   *  This function fills the 2D array \p points whose <em>(i, j)</em>th
   *  element contains the <em>i</em>th coordinate of the interpolation point
   *  \f$x_j\f$. */
  virtual void
  getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const {
    throw NotImplementedError(
        "Space::getGlobalDofInterpolationPoints(): not implemented");
  }

  /** \brief Retrieve the unit vectors normal to the grid at the
   *  interpolation points of the global degrees of freedom.
   *
   *  This function fills the 2D array \p normals whose <em>(i, j)</em>th
   *  element contains the <em>i</em>th component of the unit vector normal
   *  to the grid at the interpolation point \f$x_j\f$ defined in the
   *  documentation of getGlobalDofInterpolationPoints(). */
  virtual void getNormalsAtGlobalDofInterpolationPoints(
      Matrix<CoordinateType> &normals) const {
    throw NotImplementedError(
        "Space::getNormalsAtGlobalDofInterpolationPoints(): not implemented");
  }

  /** \brief Retrieve the interpolation directions of the global degrees of
   *  freedom.
   *
   *  This function fills the 2D array \p points whose <em>(i, j)</em>th
   *  element contains the <em>i</em>th component of the vector \f$d_j\f$
   *  defined in the documentation of getGlobalDofInterpolationPoints(). */
  virtual void getGlobalDofInterpolationDirections(
      Matrix<CoordinateType> &directions) const {
    throw NotImplementedError(
        "Space::getGlobalDofInterpolationDirections(): not implemented");
  }

  /** @}
      @name DOF reference positions (functions used in ACA assembly)
      @} */

  // These functions are used only by the ACA assembler.
  // For the moment, Point will always be 3D, independently from the
  // actual dimension of the space. Once Ahmed's bemcluster is made dimension-
  // independent, we may come up with a more elegant solution.
  /** \brief Retrieve bounding boxes of global degrees of freedom.
   *
   *  \param[out] boundingBoxes
   *    Vector whose <em>i</em>th element contains the bounding box
   *    of <em>i</em>th global degree of freedom.
   *
   *  \note This function is intended as a helper for clustering algorithms
   *  used in matrix compression algorithms such as adaptive cross
   *  approximation. */
  virtual void getGlobalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const {
    throw NotImplementedError("Space::getGlobalDofBoundingBoxes(): "
                              "implementation missing");
  }

  /** \brief Retrieve bounding boxes of local degrees of freedom ordered by
   *  their flat index.
   *
   *  \param[out] boundingBoxes
   *    Vector whose <em>i</em>th element contains the bounding box
   *    the local degree of freedom with flat index <em>i</em>.
   *
   *  \note This function is intended as a helper for clustering algorithms
   *  used in matrix compression algorithms such as adaptive cross
   *  approximation. */
  virtual void getFlatLocalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const {
    throw NotImplementedError("Space::getFlatLocalDofBoundingBoxes(): "
                              "implementation missing");
  }

  /** \brief Retrieve the reference positions of global degrees of freedom.
   *
   *  \param[out] positions
   *    Vector whose <em>i</em>th element contains the coordinates
   *    of the point taken to be the "reference position" (in some sense) of
   *    <em>i</em>th global degree of freedom.
   *
   *  \note This function is intended as a helper for clustering algorithms
   *  used in matrix compression algorithms such as adaptive cross
   *  approximation. */
  virtual void getGlobalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const = 0;

  /** \brief Retrieve the reference positions of local degrees of freedom
   *  ordered by their flat index.
   *
   *  \param[out] positions
   *    Vector whose <em>i</em>th element contains the coordinates
   *    of the point taken to be the ``reference position'' (in some sense) of
   *    the local degree of freedom with flat index <em>i</em>.
   *
   *  \note This function is intended as a helper for clustering algorithms
   *  used in matrix compression algorithms such as adaptive cross
   *  approximation. */
  virtual void getFlatLocalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const = 0;

  /** \brief Retrieve the unit vectors normal to the grid at the positions of
   *  global degrees of freedom.
   *
   *  \param[out] normals
   *    Vector whose <em>i</em>th element contains the unit vector normal to
   *    the grid at the reference position of <em>i</em>th global degree of
   *    freedom. */
  virtual void
  getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
    throw NotImplementedError("Space::getGlobalDofNormals(): not implemented");
  }

  /** \brief Retrieve the unit vectors normal to the grid at the positions of
   *  local degrees of freedom.
   *
   *  \param[out] normals
   *    Vector whose <em>i</em>th element contains the unit vector normal to
   *    the grid at the reference position of the local degree of freedom with
   *    flat index <em>i</em>. */
  virtual void
  getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
    throw NotImplementedError(
        "Space::getFlatLocalDofNormals(): not implemented");
  }

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
   *    cluster to which <em>i</em>th global degree of freedom has been
   *assigned.
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
      const char *fileName,
      const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const = 0;

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
  virtual void
  dumpClusterIdsEx(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
                   DofType dofType) const;

  
  /** \brief Initialize the cluster tree associated with this space. */
  virtual void initializeClusterTree(const ParameterList& parameterList);

  /** \brief Return the cluster tree associated with this space. */
  virtual shared_ptr<const hmat::DefaultClusterTreeType> clusterTree() const;

  /** @} */
private:
  /** \cond PRIVATE */
  shared_ptr<const Grid> m_grid;
  shared_ptr<GeometryFactory> m_elementGeometryFactory;
  unsigned int m_level;
  std::unique_ptr<GridView> m_view;
  shared_ptr<const hmat::DefaultClusterTreeType> m_clusterTree;
  /** \endcond */
};

/** \relates Space
 *  \brief Get pointers to Basis objects corresponding to all elements of the
 *grid
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
void getAllShapesets(
    const Space<BasisFunctionType> &space,
    std::vector<const Fiber::Shapeset<BasisFunctionType> *> &shapesets);

/** \relates Space
 *  \brief Get pointers to Basis objects corresponding to all elements of the
 *grid
 *  on which a function space is defined.
 *
 *  \deprecated This function is deprecated. Use getAllShapesets() instead. */
template <typename BasisFunctionType>
void BEMPP_DEPRECATED
getAllBases(const Space<BasisFunctionType> &space,
            std::vector<const Fiber::Basis<BasisFunctionType> *> &bases);

/** \relates Space
 *  \brief Return the maximum polynomial order of the shapesets
 *  defined by \p space on the elements of the space's underlying grid. */
template <typename BasisFunctionType>
int maximumShapesetOrder(const Space<BasisFunctionType> &space);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType>>
constructOperatorMappingGlobalToFlatLocalDofs(
    const Space<BasisFunctionType> &space);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType>>
constructOperatorMappingFlatLocalToGlobalDofs(
    const Space<BasisFunctionType> &space);

} // namespace Bempp

#endif
