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

#ifndef piecewise_constant_scalar_space_barycentric_hpp
#define piecewise_constant_scalar_space_barycentric_hpp

#include "../common/common.hpp"

#include "../grid/grid_view.hpp"
#include "scalar_space.hpp"
#include "../common/types.hpp"
#include "../fiber/constant_scalar_shapeset.hpp"
#include "../grid/grid_segment.hpp"
#include "../common/shared_ptr.hpp"

#include <tbb/mutex.h>

#include <map>
#include <memory>

namespace Bempp {

/** \cond FORWARD_DECL */
class GridView;
/** \endcond */

/** \ingroup space
 *  \brief Space of piecewise constant scalar functions. */
template <typename BasisFunctionType>
class PiecewiseConstantScalarSpaceBarycentric
    : public ScalarSpace<BasisFunctionType> {
public:
  typedef
      typename ScalarSpace<BasisFunctionType>::CoordinateType CoordinateType;

  /** \brief Constructor.
   *
   *  Construct a space of piecewise constant scalar functions
   *  defined on the barycentric refinement of the grid \p grid.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  explicit PiecewiseConstantScalarSpaceBarycentric(
      const shared_ptr<const Grid> &grid);

  /** \brief Constructor.
   *
   *  Construct a space of piecewise constant scalar functions defined on the
   *  elements of the grid \p grid belonging to the segment \p segment.
   *
   *  An exception is thrown if \p grid is a null pointer.
   */
  PiecewiseConstantScalarSpaceBarycentric(
      const shared_ptr<const Grid> &grid, const GridSegment &segment,
      bool strictlyOnSegment = false,
      bool putDofsOnBoundaries = true);

  virtual shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;
  virtual bool isDiscontinuous() const;

  virtual bool isBarycentric() const { return true; }

  shared_ptr<const Space<BasisFunctionType>> barycentricSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;

  virtual int domainDimension() const;
  virtual int codomainDimension() const;

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

  virtual bool spaceIsCompatible(const Space<BasisFunctionType> &other) const;

  virtual SpaceIdentifier spaceIdentifier() const {
    return PIECEWISE_CONSTANT_SCALAR_BARYCENTRIC;
  }

  virtual size_t globalDofCount() const;
  virtual size_t flatLocalDofCount() const;
  virtual void getGlobalDofs(const Entity<0> &element,
                             std::vector<GlobalDofIndex> &dofs) const;
  virtual void
  global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                   std::vector<std::vector<LocalDof>> &localDofs) const;
  virtual void
  flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &globalDofs,
                      std::vector<LocalDof> &localDofs) const;

  virtual void
  getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const;
  virtual void getNormalsAtGlobalDofInterpolationPoints(
      Matrix<CoordinateType> &normals) const;

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
  Fiber::ConstantScalarShapeset<BasisFunctionType> m_shapeset;
  std::vector<std::vector<GlobalDofIndex>> m_local2globalDofs;
  std::vector<std::vector<LocalDof>> m_global2localDofs;
  std::vector<LocalDof> m_flatLocal2localDofs;
  GridSegment m_segment;
  bool m_strictlyOnSegment;
  bool m_putDofsOnBoundaries;
  std::unique_ptr<GridView> m_view;
  shared_ptr<const Grid> m_originalGrid;
  mutable shared_ptr<Space<BasisFunctionType>> m_discontinuousSpace;
  mutable tbb::mutex m_discontinuousSpaceMutex;
  mutable Matrix<int> m_sonMap;
};

/** \brief Define a PiecewiseConstantScalarSpaceBarycentric that has an
 * update method for grid refinement. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid);

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid,
    const std::vector<int> &domains);

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid,
    int domain);

} // namespace Bempp

#endif
