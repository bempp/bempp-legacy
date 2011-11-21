// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_grid_decl_hpp
#define bempp_grid_decl_hpp

#include "id_set_decl.hpp"

#include <memory>

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;
class GridView;

/** \brief Abstract wrapper of a grid.

    Functions related to parallelisation are not wrapped yet.
 */
class Grid
{
public:
    /** Destructor */
    virtual ~Grid() {
    }

    /** @name Views @{ */

    /** \brief View for a grid level for All_Partition. */
    virtual std::auto_ptr<GridView> levelView(int level) const = 0;

    /** \brief View for the leaf grid for All_Partition. */
    virtual std::auto_ptr<GridView> leafView() const = 0;

    /** @} @name Id sets @{ */

    /// Reference to the grid's global id set
    virtual const IdSet& globalIdSet() const = 0;

    /** @} @name Adaptivity and grid refinement @{ */

    /** \brief Refine the grid refCount times using the default
     refinement rule.

     This behaves like marking all elements for refinement and then
     calling preAdapt, adapt and postAdapt. The state after
     globalRefine is comparable to the state after postAdapt.
     */
    virtual void globalRefine(int refCount) = 0;

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

     \param[in] refCount Number of subdivisions that should be
     applied. Negative value means coarsening.

     \param[in] e        Entity that should be marked

     \return true if Entity was marked, false otherwise.
     */
    virtual bool mark(int refCount, const Entity<0>& e) = 0;

    /** \brief returns adaptation mark for given entity

     \param[in] e   Entity for which adaptation mark should be determined

     \return int adaptation mark currently set for given Entity e
     */
    virtual int getMark(const Entity<0>& e) const = 0;

    /** \brief To be called after entities have been marked and before
     adapt() is called.

     This sets the mightVanish flags of the elements for the next adapt call.

     \return true if an entity may be coarsened during a subsequent
     adapt(), false otherwise.
     */
    virtual bool preAdapt() = 0;

    /** \brief Refine all positive marked leaf entities, coarsen all
     negative marked entities if possible

     \return true if a least one entity was refined

     The complete adaptation process works as follows:

     - mark entities with the mark() method
     - call preAdapt()
     - if preAdapt() returned true: possibly save current solution
     - call adapt()
     - if adapt() returned true: possibly interpolate the (saved) solution
     - call postAdapt()
     */
    virtual bool adapt() = 0;

    /** \brief To be called after grid has been adapted and
     information left over by the adaptation has been processed.

     This removes the isNew flags of the elements from the last
     adapt call.
     */
    virtual void postAdapt() = 0;

    /** @} @name Others @{ */

    /** Dimension of the space containing the grid */
    virtual int dimWorld() const = 0;

    /** \brief Maximum level defined in this grid.

     Levels are numbered 0 ... maxLevel with 0 the coarsest level.
     */
    virtual int maxLevel() const = 0;

    /** Number of boundary segments within the macro grid */
    virtual size_t numBoundarySegments() const = 0;

    /** @} */
};

/**
 Wrapper of a Dune surface grid of type DuneGrid.

 \internal The wrapper holds a pointer to a Dune Grid object. The
 member variable m_owns_dune_grid determines whether this object is
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
    /** Underlying Dune grid's type*/
    typedef DuneGrid DuneGridType;

    /** Wrap a new Dune grid object (deleted in the destructor) */
    ConcreteGrid() :
        m_dune_grid(new DuneGrid), m_owns_dune_grid(true), m_global_id_set(
            &m_dune_grid->globalIdSet()) {
    }

    /** Wrap an existing Dune grid object
     @param[in] own delete dune_grid in the destructor?

     \todo Perhaps this constructor should be removed; the user would then be
     forced to work completely in the framework of Bempp (all Dune grids managed
     by Bempp would be created by Bempp itself). Note that all other Dune objects
     created by this grid are automatically encapsulated by Bempp.
     */

    explicit ConcreteGrid(DuneGrid* dune_grid, bool own = false) :
        m_dune_grid(dune_grid), m_owns_dune_grid(own), m_global_id_set(
            dune_grid ? &dune_grid->globalIdSet() : 0) { // safety net
    }

    /// Destructor. Deletes the internal Dune Grid object.
    ~ConcreteGrid() {
        if (m_owns_dune_grid)
            delete m_dune_grid;
    }

    /// Access the underlying Dune grid object. Use at your own risk!
    const DuneGrid& duneGrid() const {
        return *m_dune_grid;
    }
    DuneGrid& duneGrid() {
        return *m_dune_grid;
    }

    virtual std::auto_ptr<GridView> levelView(int level) const;

    virtual std::auto_ptr<GridView> leafView() const;

    virtual const IdSet& globalIdSet() const {
        return m_global_id_set;
    }

    virtual void globalRefine(int refCount) {
        m_dune_grid->globalRefine(refCount);
    }

    virtual bool mark(int refCount, const Entity<0>& e);

    virtual int getMark(const Entity<0>& e) const;

    virtual bool preAdapt() {
        return m_dune_grid->preAdapt();
    }

    virtual bool adapt() {
        return m_dune_grid->adapt();
    }

    virtual void postAdapt() {
        m_dune_grid->postAdapt();
    }

    virtual int dimWorld() const {
        return DuneGrid::dimensionworld;
    }

    virtual int maxLevel() const {
        return m_dune_grid->maxLevel();
    }

    virtual size_t numBoundarySegments() const {
        return m_dune_grid->numBoundarySegments();
    }

private:
    // Disable copy constructor and assignment operator
    // (unclear what to do with the pointer to the grid)
    ConcreteGrid(const ConcreteGrid&);
    ConcreteGrid& operator =(const ConcreteGrid&);
};

} // namespace Bempp

#endif
