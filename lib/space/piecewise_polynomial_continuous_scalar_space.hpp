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

#ifndef bempp_piecewise_polynomial_continuous_scalar_space_hpp
#define bempp_piecewise_polynomial_continuous_scalar_space_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"
#include "../grid/grid_segment.hpp"

#include "scalar_space.hpp"

#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp {

/** \cond FORWARD_DECL */
class GridView;
template <typename CoordinateType> struct BoundingBox;
/** \endcond */

/** \ingroup space
 *  \brief Space of continuous, piecewise polynomial scalar functions. */
template <typename BasisFunctionType>
class PiecewisePolynomialContinuousScalarSpace
    : public ScalarSpace<BasisFunctionType> {
public:
  typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
  typedef typename Space<BasisFunctionType>::ComplexType ComplexType;

  /** \brief Constructor.
   *
   *  Construct a space of continuous functions whose restrictions to
   *  elements of the grid \p grid will be polynomials of order at most \p
   *  polynomialOrder. */
  PiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid> &grid,
                                           int polynomialOrder);

  /** \brief Constructor.
   *
   *  Construct a space of continuous functions whose restrictions to
   *  elements of the grid \p grid will be polynomials of order at most \p
   *  polynomialOrder. The space will contain only the basis functions deemed
   *  to belong to the segment \p segment; specifically, vertex functions
   *  associated with vertices belonging to \p segment, edge functions
   *  associated with edges belonging to \p segment and bubble function
   *  associated with elements belonging to \p segment. If \p
   *  strictlyOnSegment is \c true, the support of the basis functions is
   *  truncated to the elements that belong to \p segment, too; in this case,
   *  the space may in fact contain discontinuous basis functions when
   *  considered on the whole \p grid, although the basis functions will be
   *  continuous when considered on the chosen grid segment.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  PiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid> &grid,
                                           int polynomialOrder,
                                           const GridSegment &segment,
                                           bool strictlyOnSegment = false);
  virtual ~PiecewisePolynomialContinuousScalarSpace();

  virtual int domainDimension() const;
  virtual int codomainDimension() const;

  virtual bool isBarycentric() const { return false; }

  virtual bool spaceIsCompatible(const Space<BasisFunctionType> &other) const;

  virtual SpaceIdentifier spaceIdentifier() const {
    return PIECEWISE_POLYNOMIAL_CONTINUOUS_SCALAR;
  }

  /** \brief Return the variant of element \p element.
   *
   *  Possible return values:
   *    - 2: one-dimensional segment,
   *    - 3: triangular element,
   *    - 4: quadrilateral element. */
  virtual ElementVariant elementVariant(const Entity<0> &element) const;
  virtual void setElementVariant(const Entity<0> &element,
                                 ElementVariant variant);

  virtual const Fiber::Shapeset<BasisFunctionType> &
  shapeset(const Entity<0> &element) const;

  virtual shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;
  virtual bool isDiscontinuous() const;

  virtual size_t globalDofCount() const;
  virtual size_t flatLocalDofCount() const;
  virtual void getGlobalDofs(const Entity<0> &element,
                             std::vector<GlobalDofIndex> &dofs) const;
  virtual void
  global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                   std::vector<std::vector<LocalDof>> &localDofs) const;
  virtual void
  flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                      std::vector<LocalDof> &localDofs) const;

  virtual void
  getGlobalDofPositions(std::vector<Point3D<CoordinateType>> &positions) const;
  virtual void getFlatLocalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const;

  virtual void getGlobalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &bboxes) const;
  virtual void getFlatLocalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &bboxes) const;

  virtual void
  getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const;
  virtual void
  getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const;

  virtual void
  dumpClusterIds(const char *fileName,
                 const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const;
  virtual void
  dumpClusterIdsEx(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
                   DofType dofType) const;

private:
  void initialize();
  void assignDofsImpl();

private:
  /** \cond PRIVATE */
  int m_polynomialOrder;
  GridSegment m_segment;
  bool m_strictlyOnSegment;
  boost::scoped_ptr<Fiber::Shapeset<BasisFunctionType>> m_triangleShapeset;
  std::unique_ptr<GridView> m_view;
  std::vector<std::vector<GlobalDofIndex>> m_local2globalDofs;
  std::vector<std::vector<LocalDof>> m_global2localDofs;
  std::vector<LocalDof> m_flatLocal2localDofs;
  size_t m_flatLocalDofCount;
  std::vector<BoundingBox<CoordinateType>> m_globalDofBoundingBoxes;
  mutable shared_ptr<Space<BasisFunctionType>> m_discontinuousSpace;
  mutable tbb::mutex m_discontinuousSpaceMutex;
  /** \endcond */
};

/** \brief Define a PiecewisePolynomialContinuousScalarSpace that has an update method for grid refinement. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid>& grid);

/** \brief Overload to define a set of domains for the space and whether the space contains boundary entities
 (\p open = true) or not. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid>& grid,
        const std::vector<int>& domains, bool open);

/** \brief Overlad. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid>& grid,
        int domain, bool open);

} // namespace Bempp

#endif
