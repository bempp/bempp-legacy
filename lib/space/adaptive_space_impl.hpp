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
#include <boost/bind.hpp>

namespace Bempp {

template <typename BasisFunctionType_, typename SpaceType>
AdaptiveSpace<BasisFunctionType_,SpaceType>::AdaptiveSpace(const shared_ptr<const Grid>& grid, const std::vector<int>& domains, bool closed):
    Space<BasisFunctionType_>(grid), 
    m_gridSegmentFactory(grid, domains, closed), m_level(0), m_grid(grid) {

        initialize();
    }

template <typename BasisFunctionType_, typename SpaceType>
AdaptiveSpace<BasisFunctionType_,SpaceType>::AdaptiveSpace(const shared_ptr<const Grid>& grid):
    Space<BasisFunctionType_>(grid), 
    m_gridSegmentFactory(grid), m_level(0), m_grid(grid) {

        initialize();
    }

template <typename BasisFunctionType_, typename SpaceType>
AdaptiveSpace<BasisFunctionType_,SpaceType>::~AdaptiveSpace()
{
    m_connection.disconnect();
}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::initialize()
{

    m_space = shared_ptr<Space<BasisFunctionType_>>(new SpaceType(m_grid, 
            m_gridSegmentFactory.update()));
    m_connection = m_grid->connect(boost::bind(&Space<BasisFunctionType>::update,this));

}

template <typename BasisFunctionType_, typename SpaceType>
AdaptiveSpace<BasisFunctionType_,SpaceType>::AdaptiveSpace(const AdaptiveSpace<BasisFunctionType_,SpaceType>& other): 
    Space<BasisFunctionType>(other), m_gridSegmentFactory(other.m_gridSegmentFactory), m_level(other.m_level), m_grid(other.m_grid), m_space(other.m_space)
{
}

template <typename BasisFunctionType_, typename SpaceType>
bool AdaptiveSpace<BasisFunctionType_,SpaceType>::isDiscontinuous() const {

    return currentSpace().isDiscontinuous();

}

template <typename BasisFunctionType_, typename SpaceType>
shared_ptr<const Space<BasisFunctionType_>> AdaptiveSpace<BasisFunctionType_,SpaceType>::discontinuousSpace(
      const shared_ptr<const Space<BasisFunctionType_>> &self) const
{

    return currentSpace().discontinuousSpace(self);

}


template <typename BasisFunctionType_, typename SpaceType>
bool AdaptiveSpace<BasisFunctionType_,SpaceType>::isBarycentric() const
{

    return currentSpace().isBarycentric();

}

template <typename BasisFunctionType_, typename SpaceType>
shared_ptr<const Grid> AdaptiveSpace<BasisFunctionType_,SpaceType>::grid() const
{

    return m_grid;

}

template <typename BasisFunctionType_, typename SpaceType>
const GridView & AdaptiveSpace<BasisFunctionType_,SpaceType>::gridView() const
{

    return currentSpace().gridView();

}

template <typename BasisFunctionType_, typename SpaceType>
int AdaptiveSpace<BasisFunctionType_,SpaceType>::domainDimension() const
{

    return currentSpace().domainDimension();

}

template <typename BasisFunctionType_, typename SpaceType>
int AdaptiveSpace<BasisFunctionType_,SpaceType>::codomainDimension() const
{

    return currentSpace().codomainDimension();

}

template <typename BasisFunctionType_, typename SpaceType>
BEMPP_DEPRECATED const Fiber::Shapeset<BasisFunctionType_> &
  AdaptiveSpace<BasisFunctionType_,SpaceType>::shapeset(const Entity<0> &element) const
{

    return currentSpace().shapeset(element);
}

template <typename BasisFunctionType_, typename SpaceType>
BEMPP_DEPRECATED const typename AdaptiveSpace<BasisFunctionType_,SpaceType>::CollectionOfBasisTransformations &
  AdaptiveSpace<BasisFunctionType_,SpaceType>::shapeFunctionValue() const
{

    return currentSpace().shapeFunctionValue();
    
}

template <typename BasisFunctionType_, typename SpaceType>
shared_ptr<const Space<BasisFunctionType_> > 
AdaptiveSpace<BasisFunctionType_,SpaceType>::barycentricSpace(
  const shared_ptr<const Space<BasisFunctionType> > &self) const
{

    return currentSpace().barycentricSpace(self);
}


template <typename BasisFunctionType_, typename SpaceType>
const typename AdaptiveSpace<BasisFunctionType_,SpaceType>::CollectionOfShapesetTransformations & 
AdaptiveSpace<BasisFunctionType_,SpaceType>::basisFunctionValue() const
{

    return currentSpace().basisFunctionValue();

}

template <typename BasisFunctionType_, typename SpaceType> 
void AdaptiveSpace<BasisFunctionType_,SpaceType>::setElementVariant(const Entity<0> &element,
                             ElementVariant variant) {
    currentSpace().setElementVariant(element, variant);

}

template <typename BasisFunctionType_, typename SpaceType>
ElementVariant AdaptiveSpace<BasisFunctionType_,SpaceType>::elementVariant(const Entity<0> &element) const
{

    return currentSpace().elementVariant(element);

}

template <typename BasisFunctionType_, typename SpaceType>
size_t AdaptiveSpace<BasisFunctionType_,SpaceType>::flatLocalDofCount() const
{

    return currentSpace().flatLocalDofCount();

}

template <typename BasisFunctionType_, typename SpaceType>
size_t AdaptiveSpace<BasisFunctionType_,SpaceType>::globalDofCount() const
{

    return currentSpace().globalDofCount();

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofs(const Entity<0> &element,
      std::vector<GlobalDofIndex> &dofs) const
{

    return currentSpace().getGlobalDofs(element,dofs);

}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofs(const Entity<0> &element, std::vector<GlobalDofIndex> &dofs, std::vector<BasisFunctionType> &localDofWeights) const 
{

    currentSpace().getGlobalDofs(element,dofs,localDofWeights);

}

template <typename BasisFunctionType_, typename SpaceType>
bool AdaptiveSpace<BasisFunctionType_,SpaceType>::gridIsIdentical(const Space<BasisFunctionType> &other) const
{

    return currentSpace().gridIsIdentical(other);

}

template <typename BasisFunctionType_, typename SpaceType>
SpaceIdentifier AdaptiveSpace<BasisFunctionType_,SpaceType>::spaceIdentifier() const
{

    return currentSpace().spaceIdentifier();

}

template <typename BasisFunctionType_, typename SpaceType>
bool
AdaptiveSpace<BasisFunctionType_,SpaceType>::spaceIsCompatible(const Space<BasisFunctionType> &other) const
{

    return currentSpace().spaceIsCompatible(other);

}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
               std::vector<std::vector<LocalDof>> &localDofs) const 
{

    currentSpace().global2localDofs(globalDofs,localDofs);

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::global2localDofs(
  const std::vector<GlobalDofIndex> &globalDofs,
  std::vector<std::vector<LocalDof>> &localDofs,
  std::vector<std::vector<BasisFunctionType>> &localDofWeights) const
{

    currentSpace().global2localDofs(globalDofs,localDofs,localDofWeights);


}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                  std::vector<LocalDof> &localDofs) const
{

    currentSpace().flatLocal2localDofs(flatLocalDofs,localDofs);


}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const
{

    currentSpace().getGlobalDofInterpolationPoints(points);

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getNormalsAtGlobalDofInterpolationPoints(
  Matrix<CoordinateType> &normals) const 
{
    currentSpace().getNormalsAtGlobalDofInterpolationPoints(normals);

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofInterpolationDirections(
  Matrix<CoordinateType> &directions) const 
{
    currentSpace().getGlobalDofInterpolationDirections(directions);
}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofBoundingBoxes(
  std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const
{
    currentSpace().getGlobalDofBoundingBoxes(boundingBoxes);
}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getFlatLocalDofBoundingBoxes(
  std::vector<BoundingBox<CoordinateType>> &boundingBoxes) const 
{
    currentSpace().getFlatLocalDofBoundingBoxes(boundingBoxes);
}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofPositions(
  std::vector<Point3D<CoordinateType>> &positions) const
{
    currentSpace().getGlobalDofPositions(positions);
}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::getFlatLocalDofPositions(
  std::vector<Point3D<CoordinateType>> &positions) const
{
    currentSpace().getFlatLocalDofPositions(positions);

}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const
{

    currentSpace().getGlobalDofNormals(normals);

}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::getFlatLocalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const
{

    currentSpace().getFlatLocalDofNormals(normals);

}

template <typename BasisFunctionType_, typename SpaceType>
BEMPP_DEPRECATED void AdaptiveSpace<BasisFunctionType_,SpaceType>::dumpClusterIds(
  const char *fileName,
  const std::vector<unsigned int> &clusterIdsOfGlobalDofs) const
{

    currentSpace().dumpClusterIds(fileName,clusterIdsOfGlobalDofs);

}

template <typename BasisFunctionType_, typename SpaceType>
void
AdaptiveSpace<BasisFunctionType_,SpaceType>::dumpClusterIdsEx(const char *fileName,
               const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
               DofType dofType) const 
{
    currentSpace().dumpClusterIdsEx(fileName,clusterIdsOfGlobalDofs,dofType);

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::initializeClusterTree(const ParameterList& parameterList) 
{

    currentSpace().initializeClusterTree(parameterList);

}

template <typename BasisFunctionType_, typename SpaceType>
void AdaptiveSpace<BasisFunctionType_,SpaceType>::update()
{

   m_space = shared_ptr<Space<BasisFunctionType_>>(new SpaceType(m_grid, m_gridSegmentFactory.update()));
   m_level = m_grid->maxLevel(); 

}

template <typename BasisFunctionType_, typename SpaceType>
int AdaptiveSpace<BasisFunctionType_,SpaceType>::currentLevel()
{

    return m_level;

}

template <typename BasisFunctionType_, typename SpaceType>
const Space<BasisFunctionType_>& AdaptiveSpace<BasisFunctionType_,SpaceType>::currentSpace() const
{

    return *m_space;

}

template <typename BasisFunctionType_, typename SpaceType>
Space<BasisFunctionType_>& AdaptiveSpace<BasisFunctionType_,SpaceType>::currentSpace()
{

   return *m_space;

}

}
