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

#include "grid.hpp"
#include "concrete_entity.hpp"
#include "concrete_geometry_factory.hpp"
#include "concrete_grid_view.hpp"
#include "concrete_id_set.hpp"

#include <memory>

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;
class GridView;

/**
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
    ConcreteIdSet<DuneGrid, typename DuneGrid::Traits::GlobalIdSet> m_global_id_set;

public:
    /** \brief Underlying Dune grid's type*/
    typedef DuneGrid DuneGridType;

    /** \brief Wrap a new Dune grid object (deleted in the destructor). */
    ConcreteGrid() :
        m_dune_grid(new DuneGrid), m_owns_dune_grid(true), m_global_id_set(
            &m_dune_grid->globalIdSet()) {
    }

    /** \brief Wrap an existing Dune grid object.

     \param[in]  dune_grid  Pointer to the Dune grid to wrap.
     \param[in]  own  If true, *dune_grid is deleted in this object's destructor.
     */
    explicit ConcreteGrid(DuneGrid* dune_grid, bool own = false) :
        m_dune_grid(dune_grid), m_owns_dune_grid(own), m_global_id_set(
            dune_grid ? &dune_grid->globalIdSet() : 0) { // safety net
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

    virtual size_t boundarySegmentCount() const {
        return m_dune_grid->numBoundarySegments();
    }

    /** @}
    @name Views
    @{ */

    virtual std::auto_ptr<GridView> levelView(int level) const {
        return std::auto_ptr<GridView>(new ConcreteGridView<typename DuneGrid::LevelGridView>(
                                           m_dune_grid->levelView(level)));
    }

    virtual std::auto_ptr<GridView> leafView() const {
        return std::auto_ptr<GridView>(new ConcreteGridView<typename DuneGrid::LeafGridView>(
                                           m_dune_grid->leafView()));
    }

    /** @}
    @name Geometry factory
    @{ */

    virtual std::auto_ptr<GeometryFactory> elementGeometryFactory() const {
        return std::auto_ptr<GeometryFactory>(
                    new ConcreteGeometryFactory<typename DuneGrid::template Codim<0>::Geometry>());
    }

    /** @}
    @name Id sets
    @{ */

    virtual const IdSet& globalIdSet() const {
        return m_global_id_set;
    }

    /** @} */
    /** @}
    @name Adaptivity and grid refinement
    @{ */

    virtual void refineGlobally(int refCount) {
        m_dune_grid->globalRefine(refCount);
    }

    virtual bool mark(int refCount, const Entity<0>& e) {
        /// FIXME: should we catch std::bad_cast or leave it to the user?
        typedef typename DuneGrid::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_grid->mark(refCount, ce.duneEntity());
    }

    virtual int getMark(const Entity<0>& e) const {
        /// FIXME: should we catch std::bad_cast or leave it to the user?
        typedef typename DuneGrid::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_grid->getMark(ce.duneEntity());
    }

    virtual bool preAdapt() {
        return m_dune_grid->preAdapt();
    }

    virtual bool adapt() {
        return m_dune_grid->adapt();
    }

    virtual void postAdapt() {
        m_dune_grid->postAdapt();
    }

private:
    // Disable copy constructor and assignment operator
    // (unclear what to do with the pointer to the grid)
    ConcreteGrid(const ConcreteGrid&);
    ConcreteGrid& operator =(const ConcreteGrid&);
};

} // namespace Bempp

#endif
