// Copyright (C) 2011-2015 by the BEM++ Authors
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

#ifndef bempp_adaptive_space_hpp
#define bempp_adaptive_space_hpp

#include "../common/common.hpp"
#include "space.hpp"
#include "../common/shared_ptr.hpp"
#include "../grid/grid_segment.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/grid.hpp"
#include <vector>

namespace Bempp {

/** \ingroup space
 *  \brief Adaptive Function space interface.
 *
 *  This class provides a wrapper to implement a hierarchy of
 *  adaptive function spaces. */
template <typename BasisFunctionType_, typename SpaceType> 
class AdaptiveSpace : public Space<BasisFunctionType_> {
public:

  typedef BasisFunctionType_ BasisFunctionType;

  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ComplexType;

  typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
      CollectionOfShapesetTransformations;

  typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
      CollectionOfBasisTransformations;

  AdaptiveSpace(const shared_ptr<const Grid>& grid, 
          const std::vector<int>& domains, bool closed);

  AdaptiveSpace(const shared_ptr<const Grid>& grid);

  AdaptiveSpace(const AdaptiveSpace& other);

  shared_ptr<const Space<BasisFunctionType>> discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const override;

  const GridView & gridView() const override;

  shared_ptr<const Grid> grid() const override;


  bool isDiscontinuous() const override;

  bool isBarycentric() const override;

  int domainDimension() const override;

  int codomainDimension() const override;

  BEMPP_DEPRECATED const Fiber::Shapeset<BasisFunctionType> &
      shapeset(const Entity<0> &element) const override;

  BEMPP_DEPRECATED const CollectionOfBasisTransformations &
      shapeFunctionValue() const override;

  shared_ptr<const Space<BasisFunctionType>> barycentricSpace(
      const shared_ptr<const Space<BasisFunctionType>> &self) const override;


  const CollectionOfShapesetTransformations & basisFunctionValue() const override;

  void setElementVariant(const Entity<0> &element,
                                 ElementVariant variant) override;

  ElementVariant elementVariant(const Entity<0> &element) const override;

  size_t flatLocalDofCount() const override;

  size_t globalDofCount() const override;

  void getGlobalDofs(const Entity<0> &element,
          std::vector<GlobalDofIndex> &dofs) const override;

  void
  getGlobalDofs(const Entity<0> &element, std::vector<GlobalDofIndex> &dofs, std::vector<BasisFunctionType> &localDofWeights) const override;

  bool gridIsIdentical(const Space<BasisFunctionType> &other) const override;

  SpaceIdentifier spaceIdentifier() const override;
  
  bool
  spaceIsCompatible(const Space<BasisFunctionType> &other) const override;

  void
  global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                   std::vector<std::vector<LocalDof>> &localDofs) const override;

  void global2localDofs(
      const std::vector<GlobalDofIndex> &globalDofs,
      std::vector<std::vector<LocalDof>> &localDofs,
      std::vector<std::vector<BasisFunctionType>> &localDofWeights) const override;

  void
  flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                      std::vector<LocalDof> &localDofs) const override;

  void
  getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const override;

  void getNormalsAtGlobalDofInterpolationPoints(
      Matrix<CoordinateType> &normals) const override;

  void getGlobalDofInterpolationDirections(
      Matrix<CoordinateType> &directions) const override;
  
  void getGlobalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const override;

  void getFlatLocalDofBoundingBoxes(
      std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const override;

  void getGlobalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const override;

  void getFlatLocalDofPositions(
      std::vector<Point3D<CoordinateType>> &positions) const override;

  void
  getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const override;

  void
  getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const override;

  BEMPP_DEPRECATED void dumpClusterIds(
      const char *fileName,
      const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const override;

  void
  dumpClusterIdsEx(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
                   DofType dofType) const override;

  void initializeClusterTree(const ParameterList& parameterList) override;

  /** \brief Adaptively refine the space. */
  void update();

  /** \brief Return the current refinement level. */
  int currentLevel();

protected:


private:
  void initialize();

  const Space<BasisFunctionType_>& currentSpace() const;
  Space<BasisFunctionType_>& currentSpace();
    

  AdaptiveGridSegmentFactory m_gridSegmentFactory;
  int m_level;
  shared_ptr<const Grid> m_grid;
  shared_ptr<Space<BasisFunctionType>> m_space;

};

}

#include "adaptive_space_impl.hpp"

#endif
