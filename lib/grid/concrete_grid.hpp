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

#ifndef bempp_concrete_grid_hpp
#define bempp_concrete_grid_hpp

#include "../common/common.hpp"
#include "../common/ensure_not_null.hpp"
#include "grid_parameters.hpp"
#include "grid_factory.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"

#include "grid.hpp"
#include "concrete_domain_index.hpp"
#include "concrete_entity.hpp"
#include "concrete_geometry_factory.hpp"
#include "concrete_grid_view.hpp"
#include "concrete_id_set.hpp"

#include <dune/grid/common/gridview.hh>

#include <memory>

namespace Bempp {

/** \ingroup grid_internal
 *\brief Permute domain indices according to the grid index ordering. */

template <typename DuneGrid>
std::vector<int>
permuteInsertionDomainIndices(const std::vector<int> &domainIndices,
                              const Dune::GridFactory<DuneGrid> &factory,
                              const DuneGrid &grid) {

  std::vector<int> output(domainIndices.size());
  auto view = grid.leafGridView();
  const auto &indexSet = view.indexSet();
  for (auto it = grid.template leafbegin<0>(); it != grid.template leafend<0>();
       ++it) {
    const typename DuneGrid::template Codim<0>::Entity &element = *it;
    int insertionIndex = factory.insertionIndex(element);
    output[indexSet.index(element)] = domainIndices[insertionIndex];
  }

  return output;
}

/** \cond FORWARD_DECL */
template <int codim> class Entity;
class GridView;
/** \endcond */

/** \ingroup grid_internal
 \brief Wrapper of a Dune surface grid of type \p DuneGrid.

 \internal The wrapper holds a pointer to a Dune Grid object. The
 member variable \p m_owns_dune_grid determines whether this object is
 deleted in destructor.
 */
template <typename DuneGrid> class ConcreteGrid : public Grid {
private:
  DuneGrid *m_dune_grid;
  bool m_owns_dune_grid;
  GridParameters::Topology m_topology;
  ConcreteIdSet<DuneGrid, typename DuneGrid::GlobalIdSet>
      m_global_id_set;
  ConcreteDomainIndex<DuneGrid> m_domain_index;
  shared_ptr<const Dune::GridFactory<DuneGrid>> m_factory;

public:
  /** \brief Underlying Dune grid's type*/
  typedef DuneGrid DuneGridType;

  /** \brief Wrap an existing Dune grid object.

   \param[in]  dune_grid  Pointer to the Dune grid to wrap.
   \param[in]  topology   The topology of the grid.
   \param[in]  domainIndices Vector of domain indices.
   \param[in]  own  If true, *dune_grid is deleted in this object's destructor.
   */
  explicit ConcreteGrid(DuneGrid *dune_grid, GridParameters::Topology topology,
                        bool own = false)
      : m_dune_grid(ensureNotNull(dune_grid)), m_owns_dune_grid(own),
        m_topology(topology), m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(
            *dune_grid,
            std::vector<int>(dune_grid->size(0 /*level*/, 0 /*codim*/),
                             0 /*default index*/)) {}

  /** \brief Wrap an existing Dune grid object.
      \param[in] dune_grid Pointer to the Dune grid to wrap.
      \param[in] topology The topology of the grid
      \param[in] own If true, *dune_grid is deleted in this object's destructor.
  */
  explicit ConcreteGrid(DuneGrid *dune_grid, GridParameters::Topology topology,
                        const std::vector<int> &domainIndices, bool own = false)
      : m_dune_grid(ensureNotNull(dune_grid)), m_owns_dune_grid(own),
        m_topology(topology), m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(*dune_grid, domainIndices) {}

  explicit ConcreteGrid(shared_ptr<Dune::GridFactory<DuneGrid>> factory,
                        GridParameters::Topology topology)
      : m_factory(factory), m_owns_dune_grid(true),
        m_dune_grid(factory->createGrid()), m_topology(topology),
        m_global_id_set(&m_dune_grid->globalIdSet()),
        m_domain_index(*m_dune_grid,
                       std::vector<int>(m_dune_grid->size(0, 0), 0)) {}

  explicit ConcreteGrid(const shared_ptr<Dune::GridFactory<DuneGrid>> &factory,
                        GridParameters::Topology topology,
                        const std::vector<int> &domainIndices)
      : m_factory(factory), m_owns_dune_grid(true),
        m_dune_grid(factory->createGrid()), m_topology(topology),
        m_global_id_set(&m_dune_grid->globalIdSet()),
        m_domain_index(*m_dune_grid,
                       permuteInsertionDomainIndices(domainIndices, *factory,
                                                     *m_dune_grid)) {}

  /** \brief Destructor. */
  ~ConcreteGrid() {
    if (m_owns_dune_grid)
      delete m_dune_grid;
  }

  /** \brief Read-only access to the underlying Dune grid object. */
  const DuneGrid &duneGrid() const { return *m_dune_grid; }

  /** \brief Return the GridFactory used for creation of the grid (if it
     exists).
             If it does not exists an empty pointer is returned. */

  const shared_ptr<const Dune::GridFactory<DuneGrid>> factory() const {
    return m_factory;
  }

  /** \brief Access to the underlying Dune grid object. Use at your own risk! */
  DuneGrid &duneGrid() { return *m_dune_grid; }

  /** @name Grid parameters
  @{ */

  virtual int dimWorld() const override { return DuneGrid::dimensionworld; }

  virtual int dim() const override { return DuneGrid::dimension; }

  virtual int maxLevel() const override { return m_dune_grid->maxLevel(); }

  /** @}
  @name Views
  @{ */

  virtual std::unique_ptr<GridView> levelView(size_t level) const override {
    return std::unique_ptr<GridView>(
        new ConcreteGridView<typename DuneGrid::LevelGridView>(
            m_dune_grid->levelView(level), m_domain_index));
  }

  virtual std::unique_ptr<GridView> leafView() const override {
    return std::unique_ptr<GridView>(
        new ConcreteGridView<typename DuneGrid::LeafGridView>(
            m_dune_grid->leafGridView(), m_domain_index));
  }

  /** @}
  @name Geometry factory
  @{ */

  virtual std::unique_ptr<GeometryFactory> elementGeometryFactory() const override {
    return std::unique_ptr<GeometryFactory>(new ConcreteGeometryFactory<2>());
  }

  /** @}
  @name Id sets
  @{ */

  virtual const IdSet &globalIdSet() const override { return m_global_id_set; }

  /** \brief Get the grid topology */

  virtual GridParameters::Topology topology() const override { return m_topology; }

  /** @}
  @name Refinement
  @{ */

  /** \brief Return a barycentrically refined grid based on the Leaf View and its son map */
  // 
  virtual std::pair<shared_ptr<Grid>,Matrix<int>> barycentricGridSonPair() const override {
    if (!m_barycentricGrid.get()) {
      tbb::mutex::scoped_lock lock(m_barycentricSpaceMutex);
      if (!m_barycentricGrid.get()) {

        std::unique_ptr<GridView> view = this->leafView();
        const IndexSet &index=view->indexSet();

        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;

        Matrix<double> barycentricVertices;
        Matrix<int> barycentricElementCorners;
        std::vector<int> barycentricDomainIndices;

        const size_t ent0Count=view->entityCount(0); // faces
        const size_t ent1Count=view->entityCount(1); // edges
        const size_t ent2Count=view->entityCount(2); // vertices

        barycentricVertices.conservativeResize(3,ent2Count+ent1Count+ent0Count);

        for (std::unique_ptr<EntityIterator<2>> it = view->entityIterator<2>();!it->finished();it->next()){
            const Entity<2> &entity = it->entity();
            const int ent2Number = index.entityIndex(entity);
            Matrix<double> corners;
            entity.geometry().getCorners(corners);
            for(int j=0;j!=3;++j){
                barycentricVertices(j,ent2Number) = corners(j,0);
            }
        }

        for (std::unique_ptr<EntityIterator<1>> it = view->entityIterator<1>();!it->finished();it->next()){
            const Entity<1> &entity = it->entity();
            const int ent1Number = index.entityIndex(entity);
            Matrix<double> corners;
            entity.geometry().getCorners(corners);
            for(int j=0;j!=3;++j){
                barycentricVertices(j,ent2Count+ent1Number) = (corners(j,0)+corners(j,1))/2;
            }
        }

        for (std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();!it->finished();it->next()){
            const Entity<0> &entity = it->entity();
            const int ent0Number = index.entityIndex(entity);
            Matrix<double> corners;
            entity.geometry().getCorners(corners);
            for(int j=0;j!=3;++j){
                barycentricVertices(j,ent2Count+ent1Count+ent0Number) = (corners(j,0)+corners(j,1)+corners(j,2))/3;
                }
        }


        barycentricElementCorners.conservativeResize(3,6*ent0Count);
        Matrix<int> tempSonMap;
        tempSonMap.conservativeResize(6*ent0Count,2);
        for(std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();!it->finished();it->next()){
            const Entity<0> &entity = it->entity();
            const int ent0Number = index.subEntityIndex(entity,0,0);
            barycentricElementCorners(0,6*ent0Number)=index.subEntityIndex(entity,0,2);
            barycentricElementCorners(1,6*ent0Number)=ent2Count+ent1Count+ent0Number;
            barycentricElementCorners(2,6*ent0Number)=ent2Count+index.subEntityIndex(entity,1,1);

            barycentricElementCorners(0,6*ent0Number+1)=index.subEntityIndex(entity,0,2);
            barycentricElementCorners(1,6*ent0Number+1)=ent2Count+index.subEntityIndex(entity,0,1);
            barycentricElementCorners(2,6*ent0Number+1)=ent2Count+ent1Count+ent0Number;

            barycentricElementCorners(0,6*ent0Number+2)=index.subEntityIndex(entity,1,2);
            barycentricElementCorners(1,6*ent0Number+2)=ent2Count+ent1Count+ent0Number;
            barycentricElementCorners(2,6*ent0Number+2)=ent2Count+index.subEntityIndex(entity,0,1);

            barycentricElementCorners(0,6*ent0Number+3)=index.subEntityIndex(entity,1,2);
            barycentricElementCorners(1,6*ent0Number+3)=ent2Count+index.subEntityIndex(entity,2,1);
            barycentricElementCorners(2,6*ent0Number+3)=ent2Count+ent1Count+ent0Number;

            barycentricElementCorners(0,6*ent0Number+4)=index.subEntityIndex(entity,2,2);
            barycentricElementCorners(1,6*ent0Number+4)=ent2Count+ent1Count+ent0Number;
            barycentricElementCorners(2,6*ent0Number+4)=ent2Count+index.subEntityIndex(entity,2,1);

            barycentricElementCorners(0,6*ent0Number+5)=index.subEntityIndex(entity,2,2);
            barycentricElementCorners(1,6*ent0Number+5)=ent2Count+index.subEntityIndex(entity,1,1);
            barycentricElementCorners(2,6*ent0Number+5)=ent2Count+ent1Count+ent0Number;

            tempSonMap(6*ent0Number,0) = ent0Number;
            tempSonMap(6*ent0Number+1,0) = ent0Number;
            tempSonMap(6*ent0Number+2,0) = ent0Number;
            tempSonMap(6*ent0Number+3,0) = ent0Number;
            tempSonMap(6*ent0Number+4,0) = ent0Number;
            tempSonMap(6*ent0Number+5,0) = ent0Number;
            tempSonMap(6*ent0Number,1) = 0;
            tempSonMap(6*ent0Number+1,1) = 1;
            tempSonMap(6*ent0Number+2,1) = 2;
            tempSonMap(6*ent0Number+3,1) = 3;
            tempSonMap(6*ent0Number+4,1) = 4;
            tempSonMap(6*ent0Number+5,1) = 5;
        }


        shared_ptr<Grid> newGrid = 
            GridFactory::createGridFromConnectivityArrays(
                params, barycentricVertices, barycentricElementCorners, barycentricDomainIndices);

        m_barycentricGrid = newGrid;

        m_barycentricSonMap.conservativeResize(ent0Count,6);

        std::unique_ptr<GridView> baryView = m_barycentricGrid->leafView();
        const IndexSet &baryIndex=baryView->indexSet();
        int dummy = 0;
        for (std::unique_ptr<EntityIterator<0>> it = baryView->entityIterator<0>();!it->finished();it->next()){
            const Entity<0> &entity = it->entity();
            int ent0Number = baryIndex.subEntityIndex(entity,0,0);
            int insInd = m_barycentricGrid->elementInsertionIndex(entity);
            m_barycentricSonMap(tempSonMap(insInd,0),tempSonMap(insInd,1)) = ent0Number;
        }

      }
    }
    return std::pair<shared_ptr<Grid>,Matrix<int>> (m_barycentricGrid,m_barycentricSonMap);
  }

  /** \brief Return a barycentrically refined grid based on the Leaf View */
  virtual shared_ptr<Grid> barycentricGrid() const override {
    std::pair<shared_ptr<Grid>,Matrix<int>> baryPair = this->barycentricGridSonPair();
    return std::get<0>(baryPair);
  }

  /** \brief Return the son map for the barycentrically refined grid */
  virtual Matrix<int> barycentricSonMap() const override {
    std::pair<shared_ptr<Grid>,Matrix<int>> baryPair = this->barycentricGridSonPair();
    return std::get<1>(baryPair);
  }

  /** \brief Return \p true if a barycentric refinement of this grid has
   *  been created. */
  virtual bool hasBarycentricGrid() const override {
    if (!m_barycentricGrid.get())
      return false;
    else
      return true;
  }

//  /** \brief Return \p true if this is a barycentric refinement of another grid. */
//  virtual bool isBarycentricGrid() const {
//    return isBary;
//  }

  /** \brief Get insertion index of an element. */

  unsigned int elementInsertionIndex(const Entity<0> &element) const override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;
    if (!m_factory)
      throw std::runtime_error(
          "ConcreteGrid::elementInsertionIndex():: "
          "No Grid Factory defined. Cannot get insertion index.");

    return m_factory->insertionIndex(
        static_cast<const entity_t &>(element).duneEntity());
  }

  /** \brief Get insertion index of a vertex for a 2d in 3d grid. */

  unsigned int vertexInsertionIndex(const Entity<2> &vertex) const override {
    typedef ConcreteEntity<2, typename DuneGrid::template Codim<2>::Entity>
        entity_t;
    if (!m_factory)
      throw std::runtime_error(
          "ConcreteGrid::elementInsertionIndex():: "
          "No Grid Factory defined. Cannot get insertion index.");

    return m_factory->insertionIndex(
        static_cast<const entity_t &>(vertex).duneEntity());
  }

  /** \brief Pre-Adaption step */

  bool preAdapt() override {

    return m_dune_grid->preAdapt();

  }

  /** \brief Mark element for refinement. */

  bool mark(int refCount, const Entity<0>& element) override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;

    return m_dune_grid->mark(refCount, 
        static_cast<const entity_t&>(element).duneEntity());        

  }

  /** \brief Refine mesh */

  bool adapt() override {

    m_barycentricGrid.reset();
    return m_dune_grid->adapt();

  }

  /** \brief Clean up after refinement */

  void postAdapt() override {

    m_dune_grid->postAdapt();

  }

 /** \brief Refine all elements refCount times */

  void globalRefine(int refCount) override {

    m_barycentricGrid.reset();
    m_dune_grid->globalRefine(refCount);

  }

  /** \brief Return mark status of element. */
  
  int getMark(const Entity<0>& element) const override {
    typedef ConcreteEntity<0, typename DuneGrid::template Codim<0>::Entity>
        entity_t;

    return m_dune_grid->getMark(static_cast<const entity_t&>(element).duneEntity());
  }

//  /** \brief set father of barycentric refinement */
//  virtual void setBarycentricFather(shared_ptr<Grid> fatherGrid){
//    isBary=true;
//    m_barycentricFatherGrid = fatherGrid;
//  }

//  /** \brief get father of barycentric refinement */
//  virtual shared_ptr<Grid> getBarycentricFather(){
//    if(this->isBarycentricGrid()) {return m_barycentricFatherGrid;}
//    else{throw std::runtime_error("Grid is not a barycentric grid.");}
//  }


  /** @}
   */
private:
  // Disable copy constructor and assignment operator
  // (unclear what to do with the pointer to the grid)
  ConcreteGrid(const ConcreteGrid &);
  ConcreteGrid &operator=(const ConcreteGrid &);
  mutable Matrix<int> m_barycentricSonMap;
  mutable shared_ptr<Grid> m_barycentricGrid;
  mutable tbb::mutex m_barycentricSpaceMutex;
//  bool isBary=false;
//  mutable <shared_ptr<Grid> m_barycentricFatherGrid;
};

} // namespace Bempp

#endif
