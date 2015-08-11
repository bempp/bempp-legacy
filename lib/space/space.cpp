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

#include "space.hpp"
#include "../common/eigen_support.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/basis.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/geometry_factory.hpp"

#include "../hmat/geometry_data_type.hpp"
#include "../hmat/geometry.hpp"
#include "../assembly/hmat_interface.hpp"

namespace Bempp {

namespace {

template <typename BasisFunctionType>
void constructGlobalToFlatLocalDofsMappingVectors(
    const Space<BasisFunctionType> &space, std::vector<int> &rows,
    std::vector<int> &cols, std::vector<double> &values) {
  const int ldofCount = space.flatLocalDofCount();

  const GridView &view = space.gridView();
  const IndexSet &indexSet = view.indexSet();
  const size_t elementCount = view.entityCount(0);

  std::vector<std::vector<GlobalDofIndex>> gdofs(elementCount);
  std::vector<std::vector<BasisFunctionType>> ldofWeights(elementCount);
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    space.getGlobalDofs(e, gdofs[index], ldofWeights[index]);
    it->next();
  }

  rows.clear();
  cols.clear();
  rows.reserve(ldofCount);
  cols.reserve(ldofCount);

  size_t flatLdofIndex = 0;
  for (size_t e = 0; e < gdofs.size(); ++e) {
    for (size_t v = 0; v < gdofs[e].size(); ++v) {
      int gdofIndex = gdofs[e][v];
      if (gdofIndex >= 0) {
        rows.push_back(flatLdofIndex);
        cols.push_back(gdofIndex);
        ++flatLdofIndex;
      }
    }
  }
  if (rows.size() != ldofCount || cols.size() != ldofCount)
    throw std::runtime_error(
        "constructGlobalToFlatLocalDofsMappingVectors(): "
        "internal error: the number of local DOFs is different from "
        "expected. Report this problem to BEM++ developers");

  std::vector<double> tmp(ldofCount, 1.);
  values.swap(tmp);
}

template <typename BasisFunctionType>
shared_ptr<RealSparseMatrix> constructGlobalToFlatLocalDofsMappingSparseMatrix(
    const Space<BasisFunctionType> &space) {
  std::vector<int> rows, cols;
  std::vector<double> values;
  constructGlobalToFlatLocalDofsMappingVectors(space, rows, cols, values);

  assert(rows.size() == cols.size());
  assert(cols.size() == values.size());
  const size_t entryCount = rows.size();

  const size_t rowCount = space.flatLocalDofCount();
  const size_t columnCount = space.globalDofCount();

  shared_ptr<RealSparseMatrix> result =
      boost::make_shared<RealSparseMatrix>(rowCount, columnCount);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(entryCount);

  for (size_t i = 0; i < entryCount; ++i)
    triplets.push_back(Eigen::Triplet<double>(rows[i], cols[i], values[i]));
  result->setFromTriplets(triplets.begin(), triplets.end());

  return result;
}

} // namespace

template <typename BasisFunctionType>
Space<BasisFunctionType>::Space(const shared_ptr<const Grid> &grid)
    : m_grid(grid),
      m_elementGeometryFactory(grid->elementGeometryFactory().release()),
      m_view(grid->leafView()) {
  if (!grid)
    throw std::invalid_argument("Space::Space(): grid must not be a null "
                                "pointer");
}

template <typename BasisFunctionType>
Space<BasisFunctionType>::Space(const Space<BasisFunctionType> &other)
    : m_grid(other.m_grid),
      m_elementGeometryFactory(other.m_elementGeometryFactory),
      m_view(other.m_grid->levelView(other.m_level)) {}
template <typename BasisFunctionType> Space<BasisFunctionType>::~Space() {}

template <typename BasisFunctionType>
Space<BasisFunctionType> &Space<BasisFunctionType>::
operator=(const Space<BasisFunctionType> &other) {
  m_grid = other.m_grid;
  m_view = m_grid->levelView(m_level);
  m_elementGeometryFactory = other.m_elementGeometryFactory;
  return *this;
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::assignDofs() {}

template <typename BasisFunctionType>
bool Space<BasisFunctionType>::dofsAssigned() const {
  return true;
}

template <typename BasisFunctionType>
int Space<BasisFunctionType>::gridDimension() const {
  return m_grid->dim();
}

template <typename BasisFunctionType>
int Space<BasisFunctionType>::worldDimension() const {
  return m_grid->dimWorld();
}

template <typename BasisFunctionType>
const GridView &Space<BasisFunctionType>::gridView() const {
  return *m_view;
}

template <typename BasisFunctionType>
bool Space<BasisFunctionType>::gridIsIdentical(
    const Space<BasisFunctionType> &other) const {
  if (this->isBarycentric() != other.isBarycentric()) {
    return false;
  } else {
    return m_grid == other.m_grid;
  }
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
Space<BasisFunctionType>::barycentricSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
  throw std::runtime_error(
      "Space::barycentricSpace():"
      "This method is not implemented for this Space type.");
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs) const {
  throw std::runtime_error(
      "Space::getGlobalDofs(): the variant of this function "
      "with two arguments is deprecated and not "
      "implemented in new Space subclasses. Use the "
      "three-argument variant instead");
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::update()
{
    throw std::runtime_error("Space::update(): Not implemented.");

}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs,
    std::vector<BasisFunctionType> &localDofWeights) const {
  getGlobalDofs(element, dofs);
  localDofWeights.resize(dofs.size(), 1.);
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::global2localDofs(
    const std::vector<GlobalDofIndex> &globalDofs,
    std::vector<std::vector<LocalDof>> &localDofs) const {
  throw std::runtime_error(
      "Space::global2localDofs(): the variant of this function "
      "with two arguments is deprecated and not "
      "implemented in new Space subclasses. Use the "
      "three-argument variant instead");
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::global2localDofs(
    const std::vector<GlobalDofIndex> &globalDofs,
    std::vector<std::vector<LocalDof>> &localDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) const {
  global2localDofs(globalDofs, localDofs);
  localDofWeights.resize(localDofs.size());
  for (size_t igdof = 0; igdof < localDofs.size(); ++igdof)
    localDofWeights[igdof].resize(localDofs[igdof].size(),
                                  static_cast<BasisFunctionType>(1.));
}

template <typename BasisFunctionType>
boost::signals2::connection Space<BasisFunctionType>::connect(const std::function<void()>& f) const
{

    return m_spaceUpdateSignal.connect(f);

}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::sendUpdateSignal() const
{

    m_spaceUpdateSignal();

}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::initializeClusterTree(
        const ParameterList& parameterList)
{
   if (m_clusterTree)
      throw std::runtime_error("Space::initializeClusterTree(): "
                               "ClusterTree is already initialized.");

   hmat::Geometry geometry;

   auto interface = shared_ptr<hmat::GeometryInterface>(
           new SpaceHMatGeometryInterface<BasisFunctionType>(*this));
   hmat::fillGeometry(geometry,*interface);


   auto minBlockSize =
       parameterList.template get<int>("options.hmat.minBlockSize");
    m_clusterTree.reset(new hmat::DefaultClusterTreeType(
                geometry,minBlockSize));
    

}

template <typename BasisFunctionType>
shared_ptr<const hmat::DefaultClusterTreeType> 
Space<BasisFunctionType>::clusterTree() const{

    return m_clusterTree;
}

template <typename BasisFunctionType>
void getAllShapesets(
    const Space<BasisFunctionType> &space,
    std::vector<const Fiber::Shapeset<BasisFunctionType> *> &shapesets) {
  std::unique_ptr<GridView> view = space.grid()->leafView();
  const Mapper &mapper = view->elementMapper();
  const int elementCount = view->entityCount(0);

  shapesets.resize(elementCount);

  std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    shapesets[mapper.entityIndex(e)] = &space.shapeset(e);
    it->next();
  }
}

template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType> &space,
                 std::vector<const Fiber::Basis<BasisFunctionType> *> &bases) {
  std::vector<const Fiber::Shapeset<BasisFunctionType> *> shapesets;
  getAllShapesets(space, shapesets);
  bases.resize(shapesets.size());
  for (size_t i = 0; i < bases.size(); ++i) {
    bases[i] =
        dynamic_cast<const Fiber::Basis<BasisFunctionType> *>(shapesets[i]);
    if (!bases[i])
      throw std::runtime_error("getAllBases(): not all Shapeset objects "
                               "could be cast to Basis objects");
  }
}

template <typename BasisFunctionType>
int maximumShapesetOrder(const Space<BasisFunctionType> &space) {
  const GridView &view = space.gridView();
  int maxOrder = 0;
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();

  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    maxOrder = std::max(maxOrder, space.shapeset(e).order());
    it->next();
  }
  return maxOrder;
}


template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType>>
constructOperatorMappingGlobalToFlatLocalDofs(
    const Space<BasisFunctionType> &space) {
  shared_ptr<RealSparseMatrix> mat =
      constructGlobalToFlatLocalDofsMappingSparseMatrix(space);
  return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType>>(mat);
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType>>
constructOperatorMappingFlatLocalToGlobalDofs(
    const Space<BasisFunctionType> &space) {
  shared_ptr<RealSparseMatrix> mat =
      constructGlobalToFlatLocalDofsMappingSparseMatrix(space);
  return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType>>(
      mat, NO_SYMMETRY, TRANSPOSE);
}


BEMPP_GCC_DIAG_OFF(deprecated - declarations);

template <typename BasisFunctionType>
void Space<BasisFunctionType>::dumpClusterIdsEx(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
    DofType dofType) const {
  if (dofType == GLOBAL_DOFS)
    return dumpClusterIds(fileName, clusterIdsOfGlobalDofs);
  else if (dofType == FLAT_LOCAL_DOFS)
    throw std::runtime_error("Space::dumpClusterIdsEx(): "
                             "dumping of flat local DOF not supported");
  else
    throw std::invalid_argument("Space::dumpClusterIdsEx(): "
                                "invalid DOF type");
}

BEMPP_GCC_DIAG_ON(deprecated - declarations);

#define INSTANTIATE_FUNCTIONS_DEPENDENT_ON_BASIS(BASIS)                        \
  template void getAllBases(const Space<BASIS> &space,                         \
                            std::vector<const Fiber::Basis<BASIS> *> &bases);  \
  template void getAllShapesets(                                               \
      const Space<BASIS> &space,                                               \
      std::vector<const Fiber::Shapeset<BASIS> *> &bases);                     \
  template int maximumShapesetOrder(const Space<BASIS> &space)

#define INSTANTIATE_constructOperators(BASIS, RESULT)                          \
  template shared_ptr<DiscreteSparseBoundaryOperator<RESULT>>                  \
  constructOperatorMappingGlobalToFlatLocalDofs(const Space<BASIS> &space);    \
  template shared_ptr<DiscreteSparseBoundaryOperator<RESULT>>                  \
  constructOperatorMappingFlatLocalToGlobalDofs(const Space<BASIS> &space)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Space);

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FUNCTIONS_DEPENDENT_ON_BASIS);
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_constructOperators);

} // namespace Bempp
