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

#include "adaptive_space.hpp"
#include "../grid/grid.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/grid_view.hpp"

namespace Bempp {

template <typename BasisFunctionType_>
AdaptiveSpace<BasisFunctionType_>::AdaptiveSpace(const shared_ptr<const Grid>& grid, const std::vector<int>& domains, bool closed):
    Space<BasisFunctionType_>(grid), 
    m_gridSegmentFactory(grid, domains, closed), m_level(0), m_grid(grid) {
    }

template <typename BasisFunctionType_>
AdaptiveSpace<BasisFunctionType_>::AdaptiveSpace(const shared_ptr<const Grid>& grid):
    Space<BasisFunctionType_>(grid), 
    m_gridSegmentFactory(grid), m_level(0), m_grid(grid) {
    }

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::initialize()
{

    m_space = createNewSpace(m_grid, 
            m_gridSegmentFactory.update());

}

template <typename BasisFunctionType_>
AdaptiveSpace<BasisFunctionType_>::AdaptiveSpace(const AdaptiveSpace<BasisFunctionType_>& other): 
    Space<BasisFunctionType>(other), m_gridSegmentFactory(other.m_gridSegmentFactory), m_level(other.m_level), m_grid(other.m_grid), m_space(other.m_space)
{
}

template <typename BasisFunctionType_>
bool AdaptiveSpace<BasisFunctionType_>::isDiscontinuous() const {

    return currentSpace().isDiscontinuous();

}

template <typename BasisFunctionType_>
shared_ptr<const Space<BasisFunctionType_>> AdaptiveSpace<BasisFunctionType_>::discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType_>> &self) const
{

    return currentSpace().discontinuousSpace(self);

}


template <typename BasisFunctionType_>
bool AdaptiveSpace<BasisFunctionType_>::isBarycentric() const
{

    return currentSpace().isBarycentric();

}

template <typename BasisFunctionType_>
shared_ptr<const Grid> AdaptiveSpace<BasisFunctionType_>::grid() const
{

    return m_grid;

}

template <typename BasisFunctionType_>
const GridView & AdaptiveSpace<BasisFunctionType_>::gridView() const
{

    return currentSpace().gridView();

}

template <typename BasisFunctionType_>
int AdaptiveSpace<BasisFunctionType_>::domainDimension() const
{

    return currentSpace().domainDimension();

}

template <typename BasisFunctionType_>
int AdaptiveSpace<BasisFunctionType_>::codomainDimension() const
{

    return currentSpace().codomainDimension();

}

template <typename BasisFunctionType_>
BEMPP_DEPRECATED const Fiber::Shapeset<BasisFunctionType_> &
  AdaptiveSpace<BasisFunctionType_>::shapeset(const Entity<0> &element) const
{

    return currentSpace().shapeset(element);
}

template <typename BasisFunctionType_>
BEMPP_DEPRECATED const typename AdaptiveSpace<BasisFunctionType_>::CollectionOfBasisTransformations &
  AdaptiveSpace<BasisFunctionType_>::shapeFunctionValue() const
{

    return currentSpace().shapeFunctionValue();
    
}

template <typename BasisFunctionType_>
shared_ptr<const Space<BasisFunctionType_> > 
AdaptiveSpace<BasisFunctionType_>::barycentricSpace(
  const shared_ptr<const Space<BasisFunctionType> > &self) const
{

    return currentSpace().barycentricSpace(self);
}


template <typename BasisFunctionType_>
const typename AdaptiveSpace<BasisFunctionType_>::CollectionOfShapesetTransformations & 
AdaptiveSpace<BasisFunctionType_>::basisFunctionValue() const
{

    return currentSpace().basisFunctionValue();

}

template <typename BasisFunctionType_> 
void AdaptiveSpace<BasisFunctionType_>::setElementVariant(const Entity<0> &element,
                             ElementVariant variant) {
    currentSpace().setElementVariant(element, variant);

}

template <typename BasisFunctionType_>
ElementVariant AdaptiveSpace<BasisFunctionType_>::elementVariant(const Entity<0> &element) const
{

    return currentSpace().elementVariant(element);

}

template <typename BasisFunctionType_>
size_t AdaptiveSpace<BasisFunctionType_>::flatLocalDofCount() const
{

    return currentSpace().flatLocalDofCount();

}

template <typename BasisFunctionType_>
size_t AdaptiveSpace<BasisFunctionType_>::globalDofCount() const
{

    return currentSpace().globalDofCount();

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getGlobalDofs(const Entity<0> &element,
      std::vector<GlobalDofIndex> &dofs) const
{

    return currentSpace().getGlobalDofs(element,dofs);

}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::getGlobalDofs(const Entity<0> &element, std::vector<GlobalDofIndex> &dofs, std::vector<BasisFunctionType> &localDofWeights) const 
{

    currentSpace().getGlobalDofs(element,dofs,localDofWeights);

}

template <typename BasisFunctionType_>
bool AdaptiveSpace<BasisFunctionType_>::gridIsIdentical(const Space<BasisFunctionType> &other) const
{

    return currentSpace().gridIsIdentical(other);

}

template <typename BasisFunctionType_>
SpaceIdentifier AdaptiveSpace<BasisFunctionType_>::spaceIdentifier() const
{

    return currentSpace().spaceIdentifier();

}

template <typename BasisFunctionType_>
bool
AdaptiveSpace<BasisFunctionType_>::spaceIsCompatible(const Space<BasisFunctionType> &other) const
{

    return currentSpace().spaceIsCompatible(other);

}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
               std::vector<std::vector<LocalDof>> &localDofs) const 
{

    currentSpace().global2localDofs(globalDofs,localDofs);

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::global2localDofs(
  const std::vector<GlobalDofIndex> &globalDofs,
  std::vector<std::vector<LocalDof>> &localDofs,
  std::vector<std::vector<BasisFunctionType>> &localDofWeights) const
{

    currentSpace().global2localDofs(globalDofs,localDofs,localDofWeights);


}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                  std::vector<LocalDof> &localDofs) const
{

    currentSpace().flatLocal2localDofs(flatLocalDofs,localDofs);


}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const
{

    currentSpace().getGlobalDofInterpolationPoints(points);

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getNormalsAtGlobalDofInterpolationPoints(
  Matrix<CoordinateType> &normals) const 
{
    currentSpace().getNormalsAtGlobalDofInterpolationPoints(normals);

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getGlobalDofInterpolationDirections(
  Matrix<CoordinateType> &directions) const 
{
    currentSpace().getGlobalDofInterpolationDirections(directions);
}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getGlobalDofBoundingBoxes(
  std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const
{
    currentSpace().getGlobalDofBoundingBoxes(boundingBoxes);
}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getFlatLocalDofBoundingBoxes(
  std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const 
{
    currentSpace().getFlatLocalDofBoundingBoxes(boundingBoxes);
}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getGlobalDofPositions(
  std::vector<Point3D<CoordinateType>> &positions) const
{
    currentSpace().getGlobalDofPositions(positions);
}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::getFlatLocalDofPositions(
  std::vector<Point3D<CoordinateType>> &positions) const
{
    currentSpace().getFlatLocalDofPositions(positions);

}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const
{

    currentSpace().getGlobalDofNormals(normals);

}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const
{

    currentSpace().getFlatLocalDofNormals(normals);

}

template <typename BasisFunctionType_>
BEMPP_DEPRECATED void AdaptiveSpace<BasisFunctionType_>::dumpClusterIds(
  const char *fileName,
  const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const
{

    currentSpace().dumpClusterIds(fileName,clusterIdsOfGlobalDofs);

}

template <typename BasisFunctionType_>
void
AdaptiveSpace<BasisFunctionType_>::dumpClusterIdsEx(const char *fileName,
               const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
               DofType dofType) const 
{
    currentSpace().dumpClusterIdsEx(fileName,clusterIdsOfGlobalDofs,dofType);

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::initializeClusterTree(const ParameterList& parameterList) 
{

    currentSpace().initializeClusterTree(parameterList);

}

template <typename BasisFunctionType_>
void AdaptiveSpace<BasisFunctionType_>::update()
{

   m_space = createNewSpace(m_grid, m_gridSegmentFactory.update());
   m_level++;

}

template <typename BasisFunctionType_>
int AdaptiveSpace<BasisFunctionType_>::currentLevel()
{

    return m_level;

}

template <typename BasisFunctionType_>
const Space<BasisFunctionType_>& AdaptiveSpace<BasisFunctionType_>::currentSpace() const
{

    return *m_space;

}

template <typename BasisFunctionType_>
Space<BasisFunctionType_>& AdaptiveSpace<BasisFunctionType_>::currentSpace()
{

   return *m_space;

}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(AdaptiveSpace);

}
