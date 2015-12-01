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

#include "piecewise_linear_discontinuous_scalar_space_barycentric.hpp"
#include "piecewise_linear_discontinuous_scalar_space.hpp"

#include "piecewise_linear_continuous_scalar_space.hpp"
#include "space_helper.hpp"

#include "../assembly/discrete_boundary_operator.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp {

namespace {

template <typename BasisFunctionType>
class LinearDiscontinuousBarycentricSpaceFactory : public SpaceFactory<BasisFunctionType> {
    public:
       shared_ptr<Space<BasisFunctionType>> create(const shared_ptr<const Grid> &grid,
                               const GridSegment &segment) const override{
           
           return shared_ptr<Space<BasisFunctionType>>(new PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>(grid, segment));
       }
           
};


}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    PiecewiseLinearDiscontinuousScalarSpaceBarycentric(
        const shared_ptr<const Grid> &grid)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(GridSegment::wholeGrid(*grid)), m_strictlyOnSegment(false),
      m_originalGrid(grid), m_linearBasisType1(Shapeset::TYPE1),
      m_linearBasisType2(Shapeset::TYPE2), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    PiecewiseLinearDiscontinuousScalarSpaceBarycentric(
        const shared_ptr<const Grid> &grid, const GridSegment &segment,
        bool strictlyOnSegment)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(segment), m_strictlyOnSegment(strictlyOnSegment),
      m_originalGrid(grid), m_linearBasisType1(Shapeset::TYPE1),
      m_linearBasisType2(Shapeset::TYPE2), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::~PiecewiseLinearDiscontinuousScalarSpaceBarycentric() {}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::initialize() {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpaceBarycentric::initialize(): "
        "only 1- and 2-dimensional grids are supported");
  m_view = this->grid()->leafView();
  assignDofsImpl();
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (!m_discontinuousSpace) {
    tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
    typedef PiecewiseLinearDiscontinuousScalarSpace<
        BasisFunctionType> DiscontinuousSpace;
    if (!m_discontinuousSpace)
      m_discontinuousSpace.reset(new DiscontinuousSpace(
          this->grid(), m_segment, m_strictlyOnSegment));
  }
  return m_discontinuousSpace;
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    barycentricSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {

  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpaceBarycentric::barycentricSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
int PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::domainDimension() const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::codomainDimension() const {
  return 1;
}

template <typename BasisFunctionType>
bool PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::isDiscontinuous() const {
  return false;
}

template <typename BasisFunctionType>
bool PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else {
    if (other.spaceIdentifier() == PIECEWISE_LINEAR_CONTINUOUS_SCALAR) {
      // Check if this grid is a barycentric representation of the other grid
      return this->grid()->isBarycentricRepresentationOf(*other.grid());
    } else {
      return false;
    }
  }
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::assignDofsImpl() {

  const int gridDim = this->domainDimension();
  std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
  const IndexSet &index = coarseView->indexSet();
  int elementCount = this->gridView().entityCount(0);
//  int vertexCount = view.entityCount(gridDim);
  int vertexCountCoarseGrid = coarseView->entityCount(gridDim);
  int elementCountCoarseGrid = coarseView->entityCount(0);


//  const GridView &view = this->gridView();
//  int elementCountCoarseGrid = viewCoarseGrid.entityCount(0);

  // Assign gdofs to grid vertices (choosing only those that belong to
  // the selected grid segment)
  std::vector<int> globalDofIndices(3 * elementCountCoarseGrid, 0);
  int globalDofCount_ = 0;
  for (int elementIndex = 0; elementIndex < 3 * elementCountCoarseGrid; ++elementIndex)
      acc(globalDofIndices, elementIndex) = globalDofCount_++;

  int flatLocalDofCount_ = 0;



  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);
  m_elementIndex2Type.resize(elementCount);

  const int element2Basis[6][3] = {
                                   {0, 1, 2},
                                   {0, 1, 2},
                                   {2, 0, 1},
                                   {2, 0, 1},
                                   {1, 2, 0},
                                   {1, 2, 0}
                                            }; // element2Basis[i][j]
                                               // is the basis fct.
                                               // associated with the
                                               // jth vertex
                                               // on element i.*/

  // Iterate over elements
  
  globalDofCount_ = 0;
  for (std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished();it->next()) {
    const Entity<0> &entity = it->entity();
    EntityIndex ent0Number = index.subEntityIndex(entity,0,0);

    // Iterate through refined elements
    for(int i=0;i!=6;++i){
      int sonIndex = m_sonMap(ent0Number,i);

      if (i % 2 == 0) {
        acc(m_elementIndex2Type, sonIndex) = Shapeset::TYPE1;
      } else {
        acc(m_elementIndex2Type, sonIndex) = Shapeset::TYPE2;
      }

      std::vector<GlobalDofIndex> &globalDofs = acc(m_local2globalDofs, sonIndex);
      globalDofs.resize(3);



      for(int j=0;j!=3;++j){
        int basisNumber = element2Basis[i][j];
        int globalDofIndex = acc(globalDofIndices, globalDofCount_ + j);
        acc(globalDofs, basisNumber) = globalDofIndex; /// THIS LINE
//        acc(globalDofs, j) = globalDofIndex; /// THIS LINE
        acc(m_global2localDofs, globalDofIndex).push_back(LocalDof(sonIndex, basisNumber));//basisNumber
        ++flatLocalDofCount_;
      }
    }
    globalDofCount_ += 3;
  }

  // Initialize the container mapping the flat local dof indices to
  // local dof indices
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  const GridView &view = this->gridView();
  const Mapper &elementMapper = view.elementMapper();
  int index = elementMapper.entityIndex(element);
  if (m_elementIndex2Type[index] == Shapeset::TYPE1)
    return m_linearBasisType1;
  else
    return m_linearBasisType2;
}

template <typename BasisFunctionType>
ElementVariant PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::elementVariant(const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.isLine())
    return 2;
  else if (type.isTriangle())
    return 3;
  else if (type.isQuadrilateral())
    return 4;
  else
    throw std::runtime_error("PiecewiseLinearScalarSpace::"
                             "elementVariant(): invalid geometry type, "
                             "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::setElementVariant(const Entity<0> &element,
                                          ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewiseLinearScalarSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::globalDofCount() const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<
    BasisFunctionType>::getGlobalDofs(const Entity<0> &element,
                                      std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i)
    localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getGlobalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofBoundingBoxes_defaultImplementation(
          *m_view, m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  // TODO: extract this loop into a private function
  const IndexSet &indexSet = m_view->indexSet();
  const int elementCount = m_view->entityCount(0);

  std::vector<Matrix<CoordinateType>> elementCorners(elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    const Geometry &geo = e.geometry();
    geo.getCorners(acc(elementCorners, index));
    it->next();
  }

  BoundingBox<CoordinateType> model;
  const CoordinateType maxCoord = std::numeric_limits<CoordinateType>::max();
  model.lbound.x = model.lbound.y = model.lbound.z = maxCoord;
  model.ubound.x = model.ubound.y = model.ubound.z = -maxCoord;

  const int flatLocalDofCount_ = m_flatLocal2localDofs.size();
  bboxes.resize(flatLocalDofCount_, model);
  for (int i = 0; i < flatLocalDofCount_; ++i) {
    const LocalDof &localDof = acc(m_flatLocal2localDofs, i);
    BoundingBox<CoordinateType> &bbox = acc(bboxes, i);
    extendBoundingBox(bbox, acc(elementCorners, localDof.entityIndex));
    setBoundingBoxReference<CoordinateType>(
        bbox, acc(elementCorners, localDof.entityIndex).col(localDof.dofIndex));
  }

#ifndef NDEBUG
  const int globalDofCount_ = globalDofCount();
  for (size_t i = 0; i < globalDofCount_; ++i) {
    assert(bboxes[i].reference.x >= bboxes[i].lbound.x);
    assert(bboxes[i].reference.y >= bboxes[i].lbound.y);
    assert(bboxes[i].reference.z >= bboxes[i].lbound.z);
    assert(bboxes[i].reference.x <= bboxes[i].ubound.x);
    assert(bboxes[i].reference.y <= bboxes[i].ubound.y);
    assert(bboxes[i].reference.z <= bboxes[i].ubound.z);
  }
#endif // NDEBUG
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = this->domainDimension();
  const int globalDofCount_ = globalDofCount();
  const int worldDim = this->grid()->dimWorld();
  normals.resize(globalDofCount_);

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = normal(dim, 0);
    it->next();
  }

  if (gridDim == 1)
    for (size_t g = 0; g < globalDofCount_; ++g) {
      normals[g].x = 0.;
      normals[g].y = 0.;
      for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
        normals[g].x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
        normals[g].y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
      }
      auto len = std::sqrt(normals[g].x * normals[g].x + normals[g].y * normals[g].y);
      normals[g].x /= len;
      normals[g].y /= len;
    }
  else // gridDim == 2
    for (size_t g = 0; g < globalDofCount_; ++g) {
      normals[g].x = 0.;
      normals[g].y = 0.;
      for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
        normals[g].x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
        normals[g].y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
        normals[g].z += elementNormals(2, m_global2localDofs[g][l].entityIndex);
      }
      auto len = std::sqrt(normals[g].x * normals[g].x + normals[g].y * normals[g].y +
              normals[g].z * normals[g].z);
      normals[g].x /= len;
      normals[g].y /= len;
      normals[g].z /= len;
    }
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = this->domainDimension();
  const int worldDim = this->grid()->dimWorld();
  normals.resize(m_flatLocal2localDofs.size());

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = normal(dim, 0);
    it->next();
  }

  if (gridDim == 1)
    for (size_t f = 0; f < m_flatLocal2localDofs.size(); ++f) {
      int elementIndex = m_flatLocal2localDofs[f].entityIndex;
      normals[f].x = elementNormals(0, elementIndex);
      normals[f].y = elementNormals(1, elementIndex);
      normals[f].z = 0.;
    }
  else // gridDim == 2
    for (size_t f = 0; f < m_flatLocal2localDofs.size(); ++f) {
      int elementIndex = m_flatLocal2localDofs[f].entityIndex;
      normals[f].x = elementNormals(0, elementIndex);
      normals[f].y = elementNormals(1, elementIndex);
      normals[f].z = elementNormals(2, elementIndex);
    }
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    dumpClusterIds(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
    dumpClusterIdsEx(const char *fileName,
                     const std::vector<unsigned int> &clusterIdsOfDofs,
                     DofType dofType) const {
  // Note: this will probably only work for spaces on full grids
  // (not on segments)

  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpaceBarycentric::"
        "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpaceBarycentric::"
        "dumpClusterIds(): incorrect dimension");

  std::unique_ptr<GridView> view = this->grid()->leafView();
  std::unique_ptr<VtkWriter> vtkWriter = view->vtkWriter();
  if (dofType == GLOBAL_DOFS) {
    Matrix<double> data(1, idCount);
    for (size_t i = 0; i < idCount; ++i)
      data(0, i) = clusterIdsOfDofs[i];
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  } else {
    Matrix<double> data(idCount, globalDofCount());
    data.fill(0.);
    size_t row = 0;
    for (size_t id = 0; id < idCount; ++id) {
      bool exists = false;
      for (size_t fldof = 0; fldof < idCount; ++fldof) {
        if (clusterIdsOfDofs[fldof] == id) {
          LocalDof ldof = m_flatLocal2localDofs[fldof];
          GlobalDofIndex gdof =
              m_local2globalDofs[ldof.entityIndex][ldof.dofIndex];
          data(row, gdof) = 1;
          exists = true;
        }
      }
      if (!exists)
        eigenRemoveRowFromMatrix(data, row);
      else
        ++row;
    }
    std::cout << "about to write" << std::endl;
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  }
}


template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewiseLinearDiscontinuousScalarSpaceBarycentric(const shared_ptr<const Grid>& grid)
{

    shared_ptr<SpaceFactory<BasisFunctionType>> factory(
            new LinearDiscontinuousBarycentricSpaceFactory<BasisFunctionType>());
    return shared_ptr<Space<BasisFunctionType>>(
            new AdaptiveSpace<BasisFunctionType>(factory, grid));

}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)   \
    template shared_ptr<Space<BASIS>> adaptivePiecewiseLinearDiscontinuousScalarSpaceBarycentric<BASIS>( \
            const shared_ptr<const Grid>&)

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseLinearDiscontinuousScalarSpaceBarycentric);

} // namespace Bempp
