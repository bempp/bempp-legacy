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

#ifndef bempp_scaled_nedelec_0_vector_space_barycentric_hpp
#define bempp_scaled_nedelec_0_vector_space_barycentric_hpp

#include "../common/common.hpp"

#include "space.hpp"

#include "dof_assignment_mode.hpp"

#include "../grid/grid_segment.hpp"
#include "../grid/grid_view.hpp"
#include "../common/types.hpp"
#include "../fiber/nedelec_0_shapeset.hpp"
#include "../fiber/nedelec_0_shapeset_barycentric.hpp"

#include <boost/scoped_ptr.hpp>
#include <map>
#include <memory>
#include <tbb/mutex.h>

namespace Bempp {

/** \cond FORWARD_DECL */
class GridView;
class Grid;
/** \endcond */

/** \ingroup space
 *  \brief Space of continuous, piecewise linear scalar functions. */
template <typename BasisFunctionType>
class ScaledNedelec0VectorSpaceBarycentric : public Space<BasisFunctionType> {
  typedef Space<BasisFunctionType> Base;

public:
  typedef typename Space<BasisFunctionType>::CoordinateType CoordinateType;
  typedef typename Space<BasisFunctionType>::ComplexType ComplexType;
  typedef typename Base::CollectionOfShapesetTransformations
      CollectionOfShapesetTransformations;
  typedef typename Base::CollectionOfBasisTransformations
      CollectionOfBasisTransformations;

  explicit ScaledNedelec0VectorSpaceBarycentric(
      const shared_ptr<const Grid> &grid, bool putDofsOnBoundaries = false);
  ScaledNedelec0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid,
                                       const GridSegment &segment,
                                       bool putDofsOnBoundaries = false,
                                       bool strictlyOnSegment=false,
                                       int dofMode = EDGE_ON_SEGMENT);
  virtual ~ScaledNedelec0VectorSpaceBarycentric();

  virtual shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const;
  virtual bool isDiscontinuous() const;

  virtual const CollectionOfShapesetTransformations &basisFunctionValue() const;

  virtual int domainDimension() const;
  virtual int codomainDimension() const;

  virtual bool isBarycentric() const { return false; }

  virtual bool spaceIsCompatible(const Space<BasisFunctionType> &other) const;

  virtual SpaceIdentifier spaceIdentifier() const {
    return SCALED_NEDELEC_0_VECTOR_BARYCENTRIC;
  }

  /** \brief Return the variant of element \p element.
   *
   *  Possible return values:
   *    - 3: triangular element,
   *    - 4: quadrilateral element. */
  virtual ElementVariant elementVariant(const Entity<0> &element) const;
  virtual void setElementVariant(const Entity<0> &element,
                                 ElementVariant variant);

  virtual const Fiber::Shapeset<BasisFunctionType> &
  shapeset(const Entity<0> &element) const;

  virtual size_t globalDofCount() const;
  virtual size_t flatLocalDofCount() const;
  virtual void getGlobalDofs(const Entity<0> &element,
                             std::vector<GlobalDofIndex> &dofs,
                             std::vector<BasisFunctionType> &dofWeights) const;
  virtual void global2localDofs(
      const std::vector<GlobalDofIndex> &globalDofs,
      std::vector<std::vector<LocalDof>> &localDofs,
      std::vector<std::vector<BasisFunctionType>> &localDofWeights) const;
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
  struct Impl;
  boost::scoped_ptr<Impl> m_impl;
  GridSegment m_segment;
  typedef Fiber::Nedelec0ShapesetBarycentric<BasisFunctionType> Shapeset;
  bool m_putDofsOnBoundaries;
  bool m_strictlyOnSegment;
  int m_dofMode;
  std::unique_ptr<GridView> m_view;
  Fiber::Nedelec0Shapeset<3, BasisFunctionType> m_triangleShapeset;
  std::vector<std::vector<GlobalDofIndex>> m_local2globalDofs;
  std::vector<std::vector<BasisFunctionType>> m_local2globalDofWeights;
  std::vector<std::vector<LocalDof>> m_global2localDofs;
  std::vector<LocalDof> m_flatLocal2localDofs;
  std::vector<BoundingBox<CoordinateType>> m_globalDofBoundingBoxes;
  mutable shared_ptr<Space<BasisFunctionType>> m_discontinuousSpace;
  mutable tbb::mutex m_discontinuousSpaceMutex;
  std::vector<typename Shapeset::BasisType> m_elementShapesets;
  mutable Matrix<int> m_sonMap;
  mutable shared_ptr<const Grid> m_originalGrid;
  Shapeset m_RTBasisType1;
  Shapeset m_RTBasisType2;

  /** \endcond */
};

/** \brief Define a ScaledNedelec0VectorSpaceBarycentric that has an update
 * method for grid refinement. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveScaledNedelec0VectorSpaceBarycentric(
    const shared_ptr<const Grid> &grid);

/** \brief Overload to define a set of domains for the space and whether the
 space contains boundary entities
 (\p open = true) or not. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveScaledNedelec0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid,
                                             const std::vector<int> &domains,
                                             bool open, bool strictly);

/** \brief Overload. */
template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveScaledNedelec0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid,
                                             int domain, bool open, bool strictly);

} // namespace Bempp

#endif
