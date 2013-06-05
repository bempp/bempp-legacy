#ifndef DUNE_FOAMGRID_LEAFITERATOR_HH
#define DUNE_FOAMGRID_LEAFITERATOR_HH

/** \file
* \brief The FoamGridLeafIterator class
*/

namespace Dune {


/** \brief Iterator over all entities of a given codimension and level of a grid.
*  \ingroup FoamGrid
*/
template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLeafIterator :
    public Dune::FoamGridEntityPointer <codim,GridImp>
{
    enum {dim      = GridImp::dimension};
    enum {dimworld = GridImp::dimensionworld};

public:

    FoamGridLeafIterator(const GridImp& grid) 
        : FoamGridEntityPointer <codim,GridImp>(nullptr),
          grid_(&grid)
    {

        /** \todo Can a make the fullRefineLevel work somehow? */
        const int fullRefineLevel = 0;
        
        const std::list<FoamGridEntityImp<dim-codim,dimworld> >& entities = Dune::get<dim-codim>(grid_->entityImps_[fullRefineLevel]);
        // The &* turns an iterator into a plain pointer
        this->virtualEntity_.setToTarget(&*entities.begin());
        levelIterator_ = entities.begin();

        if (!this->virtualEntity_.getTarget()->isLeaf())
            increment();
    }

  //! Constructor
    FoamGridLeafIterator() 
        : FoamGridEntityPointer <codim,GridImp>(nullptr),
          grid_(nullptr)
    {}

    //! prefix increment
    void increment() {
        // Increment until you find a leaf entity
        do {
            globalIncrement();

        } while (levelIterator_!=Dune::get<dim-codim>(grid_->entityImps_[grid_->maxLevel()]).end() 
                 && !this->virtualEntity_.getTarget()->isLeaf());
    }

private:

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

        // Backup current level because it may not be accessible anymore after
        // setting the pointer to the next entity.
        const int oldLevel = this->virtualEntity_.level();

        // Increment on this level
        ++levelIterator_;
        this->virtualEntity_.setToTarget(&(*levelIterator_));
        if (levelIterator_==Dune::get<dim-codim>(grid_->entityImps_[oldLevel]).end())
            this->virtualEntity_.setToTarget(nullptr);

        // If beyond the end of this level set to first of next level
        if (levelIterator_==Dune::get<dim-codim>(grid_->entityImps_[oldLevel]).end() && oldLevel < grid_->maxLevel()) {

            const std::list<FoamGridEntityImp<dim-codim,dimworld> >& entities = Dune::get<dim-codim>(grid_->entityImps_[oldLevel+1]);
            levelIterator_ = entities.begin();
            this->virtualEntity_.setToTarget(&*entities.begin());

        }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

    // This iterator derives from FoamGridEntityPointer, and that base class stores the value
    // of the iterator, i.e. the 'pointer' to the entity.  However, that pointer can not be
    // set to its successor in the level std::list, not even by magic.  Therefore we keep the
    // same information redundantly in this iterator, which can be incremented.
    typename std::list<FoamGridEntityImp<dim-codim,dimworld> >::const_iterator levelIterator_;
};


}  // namespace Dune
  
#endif
