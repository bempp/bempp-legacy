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

#include "grid.hpp"
#include "concrete_domain_index.hpp"
#include "concrete_entity.hpp"
#include "concrete_geometry_factory.hpp"
#include "concrete_grid_view.hpp"
#include "concrete_id_set.hpp"
#include <armadillo>

#include <memory>

namespace Bempp
{

/** \cond FORWARD_DECL */
template<int codim> class Entity;
class GridView;
/** \endcond */

/** \ingroup grid_internal
 \brief Wrapper of a Dune surface grid of type \p DuneGrid.

 \internal The wrapper holds a pointer to a Dune Grid object. The
 member variable \p m_owns_dune_grid determines whether this object is
 deleted in destructor.
 */
template<typename DuneGrid>
class ConcreteGrid: public Grid
{
private:
    DuneGrid* m_dune_grid;
    bool m_owns_dune_grid;
    GridParameters::Topology m_topology;
    ConcreteIdSet<DuneGrid, typename DuneGrid::Traits::GlobalIdSet> m_global_id_set;
    ConcreteDomainIndex<DuneGrid> m_domain_index;

public:
    /** \brief Underlying Dune grid's type*/
    typedef DuneGrid DuneGridType;

    /** \brief Wrap an existing Dune grid object.

     \param[in]  dune_grid  Pointer to the Dune grid to wrap.
     \param[in]  topology   The topology of the grid.
     \param[in]  domainIndices Vector of domain indices.
     \param[in]  own  If true, *dune_grid is deleted in this object's destructor.
     */
    explicit ConcreteGrid(DuneGrid* dune_grid,
                          GridParameters::Topology topology, bool own = false) :
        m_dune_grid(ensureNotNull(dune_grid)),
        m_owns_dune_grid(own),
        m_topology(topology),
        m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(*dune_grid,
                       std::vector<int>(
                           dune_grid->size(0 /*level*/, 0 /*codim*/),
                           0 /*default index*/))
    {

    }

    /** \brief Wrap an existing Dune grid object.
        \param[in] dune_grid Pointer to the Dune grid to wrap.
        \param[in] topology The topology of the grid
        \param[in] own If true, *dune_grid is deleted in this object's destructor.
    */
    explicit ConcreteGrid(DuneGrid* dune_grid,
                          GridParameters::Topology topology,
                          const std::vector<int>& domainIndices,
                          bool own = false) :
        m_dune_grid(ensureNotNull(dune_grid)),
        m_owns_dune_grid(own),
        m_topology(topology),
        m_global_id_set(&dune_grid->globalIdSet()),
        m_domain_index(*dune_grid, domainIndices)
    {
    }

    /** \brief Destructor. */
    ~ConcreteGrid() {
        if (m_owns_dune_grid)
            delete m_dune_grid;
    }

    /** \brief Read-only access to the underlying Dune grid object. */
    const DuneGrid& duneGrid() const {
        return *m_dune_grid;
    }

    /** \brief Access to the underlying Dune grid object. Use at your own risk! */
    DuneGrid& duneGrid() {
        return *m_dune_grid;
    }

    /** @name Grid parameters
    @{ */

    virtual int dimWorld() const {
        return DuneGrid::dimensionworld;
    }

    virtual int dim() const {
        return DuneGrid::dimension;
    }

    virtual int maxLevel() const {
        return m_dune_grid->maxLevel();
    }

    /** @}
    @name Views
    @{ */

    virtual std::unique_ptr<GridView> levelView(size_t level) const {
        return std::unique_ptr<GridView>(
                    new ConcreteGridView<typename DuneGrid::LevelGridView>(
                        m_dune_grid->levelView(level), m_domain_index));
    }

    virtual std::unique_ptr<GridView> leafView() const {
        return std::unique_ptr<GridView>(
                    new ConcreteGridView<typename DuneGrid::LeafGridView>(
                        m_dune_grid->leafView(), m_domain_index));
    }

    /** @}
    @name Geometry factory
    @{ */

    virtual std::unique_ptr<GeometryFactory> elementGeometryFactory() const {
        return std::unique_ptr<GeometryFactory>(
                    new ConcreteGeometryFactory<
                    typename DuneGrid::template Codim<0>::Geometry>());
    }

    /** @}
    @name Id sets
    @{ */

    virtual const IdSet& globalIdSet() const {
        return m_global_id_set;
    }

    /** \brief Get the grid topology */

    virtual GridParameters::Topology topology() const {
        return m_topology;
    }


    /** @}
    @name Refinement
    @{ */

    /** \brief Return a barycentrically refined grid based on the LeafView */
    virtual shared_ptr<Grid> barycentricGrid() const {


        if (!m_barycentricGrid.get()){
            tbb::mutex::scoped_lock lock(m_barycentricSpaceMutex);
            if (!m_barycentricGrid.get()){

                arma::Mat<double> vertices;
                arma::Mat<int> elementCorners;
                arma::Mat<char> auxData;
                std::vector<int> domainIndices;

                std::unique_ptr<GridView> view = this->leafView();

                view->getRawElementData(vertices,
                                       elementCorners,
                                       auxData,
                                       domainIndices);

                GridParameters params;
                params.topology = GridParameters::TRIANGULAR;
                shared_ptr<Grid> newGrid = GridFactory::createGridFromConnectivityArrays(params,vertices,elementCorners,
                                                                                         domainIndices);
                shared_ptr<ConcreteGrid<DuneGrid> > concreteGrid = dynamic_pointer_cast<ConcreteGrid<DuneGrid> >(newGrid);
                (concreteGrid->duneGrid()).globalBarycentricRefine(1);

                m_barycentricGrid = newGrid;

            }
        }
        return m_barycentricGrid;
    }

    /** \brief Return \p true if a barycentric refinement of this grid has
     *  been created. */
    virtual bool hasBarycentricGrid() const {
        if (!m_barycentricGrid.get())
            return false;
        else
            return true;
    }


    /** @}
     */


private:
    // Disable copy constructor and assignment operator
    // (unclear what to do with the pointer to the grid)
    ConcreteGrid(const ConcreteGrid&);
    ConcreteGrid& operator =(const ConcreteGrid&);
    mutable shared_ptr<Grid> m_barycentricGrid;
    mutable tbb::mutex m_barycentricSpaceMutex;

};

} // namespace Bempp

#endif
