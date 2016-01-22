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

#include "piecewise_linear_discontinuous_scalar_space.hpp"
//#include "piecewise_linear_discontinuous_scalar_space_barycentric.hpp"
#include "adaptive_space.hpp"

#include "space_helper.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_segment.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp {

namespace {

template <typename BasisFunctionType>
class LinearDiscontinuousSpaceFactory : public SpaceFactory<BasisFunctionType> {

    public:
       LinearDiscontinuousSpaceFactory(bool strictlyOnSegment) :
           m_strictlyOnSegment(strictlyOnSegment){}

       shared_ptr<Space<BasisFunctionType>> create(const shared_ptr<const Grid> &grid,
                               const GridSegment &segment) const override{
           
           return shared_ptr<Space<BasisFunctionType>>(new PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>(grid, segment, m_strictlyOnSegment));
       }
    private:
       bool m_strictlyOnSegment;
           
};

}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    PiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid> &grid)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid),
      m_segment(GridSegment::wholeGrid(*grid)), m_strictlyOnSegment(false) {
  initialize(m_segment, m_strictlyOnSegment);
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    PiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid> &grid,
                                            const GridSegment &segment,
                                            bool strictlyOnSegment)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid), m_segment(segment),
      m_strictlyOnSegment(strictlyOnSegment) {
  initialize(segment, strictlyOnSegment);
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<
    BasisFunctionType>::~PiecewiseLinearDiscontinuousScalarSpace() {}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::initialize(
    const GridSegment &segment, bool strictlyOnSegment) {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpace::"
        "PiecewiseLinearDiscontinuousScalarSpace(): "
        "only 1- and 2-dimensional grids are supported");
  m_view = this->grid()->leafView();
  assignDofsImpl(segment, strictlyOnSegment);
}

template <typename BasisFunctionType>
bool PiecewiseLinearDiscontinuousScalarSpace<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else
    return false;
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewiseLinearDiscontinuousScalarSpace::discontinuousSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
bool PiecewiseLinearDiscontinuousScalarSpace<
    BasisFunctionType>::isDiscontinuous() const {
  return true;
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::assignDofsImpl(
    const GridSegment &segment, bool strictlyOnSegment) {
  const int gridDim = this->domainDimension();

  const Mapper &elementMapper = m_view->elementMapper();
  const IndexSet &indexSet = m_view->indexSet();

  //    int globalDofCount_ = m_view->entityCount(this->grid()->dim());
  int elementCount = m_view->entityCount(0);

  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  // estimated number of global DOFs
  m_global2localDofs.reserve(4 * elementCount);
  // TODO: consider calling reserve(x) for each element of m_global2localDofs
  // with x being the typical number of elements adjacent to a vertex in a
  // grid of dimension gridDim

  // Iterate over elements
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  size_t globalDofCount = 0;
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    EntityIndex elementIndex = elementMapper.entityIndex(element);

    int vertexCount;
    if (gridDim == 1)
      vertexCount = element.template subEntityCount<1>();
    else // gridDim == 2
      vertexCount = element.template subEntityCount<2>();

    // List of global DOF indices corresponding to the local DOFs of the
    // current element
    std::vector<GlobalDofIndex> &globalDofs = m_local2globalDofs[elementIndex];
    globalDofs.reserve(vertexCount);
    for (int i = 0; i < vertexCount; ++i) {
      int vertexIndex = indexSet.subEntityIndex(element, i, gridDim);
      if ((strictlyOnSegment && segment.contains(0, elementIndex)) ||
          (!strictlyOnSegment && segment.contains(gridDim, vertexIndex))) {
        globalDofs.push_back(globalDofCount);
        std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, i));
        m_global2localDofs.push_back(localDofs);
        ++globalDofCount;
      } else
        globalDofs.push_back(-1);
    }
    it->next();
  }

  // Initialize the container mapping the flat local dof indices to
  // local dof indices
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      globalDofCount, m_local2globalDofs, m_flatLocal2localDofs);
}

template <typename BasisFunctionType>
size_t
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::globalDofCount()
    const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::flatLocalDofCount()
    const {
  return globalDofCount();
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i) {
    localDofs[i] = m_global2localDofs[globalDofs[i]];
    assert(localDofs[i].size() == 1);
  }
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getGlobalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  getGlobalDofPositions(positions);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofBoundingBoxes_defaultImplementation(
          *m_view, m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  getGlobalDofBoundingBoxes(bboxes);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
      *m_view, m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType>> &normals) const {
  getGlobalDofNormals(normals);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
    dumpClusterIdsEx(const char *fileName,
                     const std::vector<unsigned int> &clusterIdsOfDofs,
                     DofType dofType) const {
  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument("PiecewiseLinearDiscontinuousScalarSpace::"
                                "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument("PiecewiseLinearDiscontinuousScalarSpace::"
                                "dumpClusterIds(): incorrect dimension");

  std::unique_ptr<GridView> view = this->grid()->leafView();
  std::unique_ptr<VtkWriter> vtkWriter = view->vtkWriter();

  Matrix<double> data(idCount, globalDofCount());
  data.setZero();
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
      eigenRemoveRowFromMatrix(data, row); // very inefficient, of course
    else
      ++row;
  }
  std::cout << "about to write" << std::endl;
  vtkWriter->addVertexData(data, "ids");
  vtkWriter->write(fileName);
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid)
{

    shared_ptr<SpaceFactory<BasisFunctionType>> factory(
            new LinearDiscontinuousSpaceFactory<BasisFunctionType>(false));
    return shared_ptr<Space<BasisFunctionType>>(
            new AdaptiveSpace<BasisFunctionType>(factory, grid));

}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid,
        const std::vector<int>& domains, bool open, bool strictlyOnSegment)
{

    shared_ptr<SpaceFactory<BasisFunctionType>> factory(
            new LinearDiscontinuousSpaceFactory<BasisFunctionType>(strictlyOnSegment));
    return shared_ptr<Space<BasisFunctionType>>(
            new AdaptiveSpace<BasisFunctionType>(factory, grid, domains, open));

}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid,
        int domain, bool open, bool strictlyOnSegment)
{
    
    shared_ptr<SpaceFactory<BasisFunctionType>> factory(
            new LinearDiscontinuousSpaceFactory<BasisFunctionType>(strictlyOnSegment));
    return shared_ptr<Space<BasisFunctionType>>(
            new AdaptiveSpace<BasisFunctionType>(factory, grid,
                std::vector<int>({domain}),open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)   \
    template shared_ptr<Space<BASIS>> adaptivePiecewiseLinearDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&); \
    template shared_ptr<Space<BASIS>> adaptivePiecewiseLinearDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&, \
            const std::vector<int>&, bool, bool); \
    template shared_ptr<Space<BASIS>> adaptivePiecewiseLinearDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&, \
            int, bool, bool) 


FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseLinearDiscontinuousScalarSpace);

} // namespace Bempp
