#ifndef DUNE_FOAMGRID_HH
#define DUNE_FOAMGRID_HH

/** \file
* \brief The FoamGrid class
*/

#include <list>

#include <dune/common/collectivecommunication.hh>
#include <dune/common/tuples.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// Implementation classes
#include "foamgrid/foamgridvertex.hh"
#include "foamgrid/foamgridedge.hh"
#include "foamgrid/foamgridelements.hh"

// The components of the FoamGrid interface
#include "foamgrid/foamgridgeometry.hh"
#include "foamgrid/foamgridentity.hh"
#include "foamgrid/foamgridentitypointer.hh"
#include "foamgrid/foamgridintersectioniterators.hh"
#include "foamgrid/foamgridleveliterator.hh"
#include "foamgrid/foamgridleafiterator.hh"
#include "foamgrid/foamgridhierarchiciterator.hh"
#include "foamgrid/foamgridindexsets.hh"

namespace Dune {

// Forward declaration
template <int dimworld>
class FoamGrid;


/** \brief Encapsulates loads of types exported by FoamGrid */
template<int dimworld>
struct FoamGridFamily
{
    typedef GridTraits<
        2,   // dim
        dimworld,   // dimworld
        Dune::FoamGrid<dimworld>,
        FoamGridGeometry,
        FoamGridEntity,
        FoamGridEntityPointer,
        FoamGridLevelIterator,
        FoamGridLeafIntersection,
        FoamGridLevelIntersection,
        FoamGridLeafIntersectionIterator,
        FoamGridLevelIntersectionIterator,
        FoamGridHierarchicIterator,
        FoamGridLeafIterator,
        FoamGridLevelIndexSet< const FoamGrid<dimworld> >,
        FoamGridLeafIndexSet< const FoamGrid<dimworld> >,
        FoamGridGlobalIdSet< const FoamGrid<dimworld> >,
        unsigned int,   // global id type
        FoamGridLocalIdSet< const FoamGrid<dimworld> >,
        unsigned int,   // local id type
        CollectiveCommunication<Dune::FoamGrid<dimworld> >
            > Traits;
};



/** \brief An implementation of the Dune grid interface: a 2d simplicial grid in an n-dimensional world
 * 
* \tparam dimworld Dimension of the world space
*/
template <int dimworld>
class FoamGrid :
        public GridDefaultImplementation  <2, dimworld, double, FoamGridFamily<dimworld> >
{
    
    friend class FoamGridLevelIndexSet<const FoamGrid >;
    friend class FoamGridLeafIndexSet<const FoamGrid >;
    friend class FoamGridGlobalIdSet<const FoamGrid >;
    friend class FoamGridLocalIdSet<const FoamGrid >;
    friend class FoamGridHierarchicIterator<const FoamGrid >;
    friend class FoamGridLevelIntersectionIterator<const FoamGrid >;
    friend class FoamGridLeafIntersectionIterator<const FoamGrid >;
    friend class FoamGridLevelIntersection<const FoamGrid >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class FoamGridLeafIterator;
    
    template <class GridType_>
    friend class GridFactory;

    template<int codim_, int dim_, class GridImp_>
    friend class FoamGridEntity;

    public:
        
    /** \brief This grid is always 2-dimensional */
    enum {dimension = 2};

    //**********************************************************
    // The Interface Methods
    //**********************************************************
    
    //! type of the used GridFamily for this grid
    typedef FoamGridFamily<dimworld>  GridFamily;
    
    //! Exports various types belonging to this grid class
    typedef typename FoamGridFamily<dimworld>::Traits Traits;
    
    //! The type used to store coordinates
    typedef double ctype;
    
    /** \brief Constructor, constructs an empty grid
     */
    FoamGrid() 
    {
        std::fill(freeIdCounter_.begin(), freeIdCounter_.end(), 0);
    }
        
        //! Destructor
        ~FoamGrid()
        {
            // Delete level index sets
            for (size_t i=0; i<levelIndexSets_.size(); i++)
                if (levelIndexSets_[i])
                    delete (levelIndexSets_[i]);
        }
        
        
        //! Return maximum level defined in this grid. Levels are numbered
        //! 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const {
            return entityImps_.size()-1;;
        }
        
        
        //! Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");

            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimworld> >(Dune::get<dimension-codim>(entityImps_[level]).begin());
        }
    
        
        //! one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,All_Partition, const Dune::FoamGrid<dimworld> >(Dune::get<dimension-codim>(entityImps_[level]).end());
        }
        
        
        //! Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(Dune::GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimworld> >(Dune::get<dimension-codim>(entityImps_[level]).begin());
        }
        

        //! one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "LevelIterator in nonexisting level " << level << " requested!");
            
            return Dune::FoamGridLevelIterator<codim,PiType, const Dune::FoamGrid<dimworld> >(Dune::get<dimension-codim>(entityImps_[level]).end());
        }
        
    
        //! Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >(*this);
        }
        
    
        //! one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,All_Partition, const FoamGrid >();
        }
        
    
        //! Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >(*this);
        }
        
        
        //! one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
            return FoamGridLeafIterator<codim,PiType, const FoamGrid >();
        }
        

        /** \brief Number of grid entities per level and codim
         */
        int size (int level, int codim) const {

            // Turn dynamic index into static index
            if (codim==2)
                return Dune::get<0>(entityImps_[level]).size();
            if (codim==1)
                return Dune::get<1>(entityImps_[level]).size();
            if (codim==0)
                return Dune::get<2>(entityImps_[level]).size();

            return 0;
        }
        
        
        //! number of leaf entities per codim in this process
        int size (int codim) const{
            return leafIndexSet().size(codim);
        }
        
        
        //! number of entities per level, codim and geometry type in this process
        int size (int level, GeometryType type) const {
            return levelIndexSets_[level]->size(type);
        }
        
            
        //! number of leaf entities per codim and geometry type in this process
        int size (GeometryType type) const
        {
            return leafIndexSet().size(type);
        }
        
        /** \brief The number of boundary edges on the coarsest level */
        size_t numBoundarySegments() const
        {
            DUNE_THROW(Dune::NotImplemented, "numBoundarySegments");
            return 0;
        }
        
        /** \brief Access to the GlobalIdSet */
        const typename Traits::GlobalIdSet& globalIdSet() const{
            return globalIdSet_;
        }
        
        
        /** \brief Access to the LocalIdSet */
        const typename Traits::LocalIdSet& localIdSet() const{
            return localIdSet_;
        }
        
        
        /** \brief Access to the LevelIndexSets */
        const typename Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return *levelIndexSets_[level];
        }
        
        
        /** \brief Access to the LeafIndexSet */
        const typename Traits::LeafIndexSet& leafIndexSet() const
        {
            return leafIndexSet_;
        }
        
        
        /** @name Grid Refinement Methods */
        /*@{*/
        
        
        /** global refinement
        */
        void globalRefine (int refCount)
        {
            DUNE_THROW(NotImplemented, "globalRefine!");
        }
        
        /** \brief Mark entity for refinement
        *
        * This only works for entities of codim 0.
        * The parameter is currently ignored
        *
        * \return <ul>
        * <li> true, if marking was successful </li>
        * <li> false, if marking was not possible </li>
        * </ul>
        */
        bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
        {
            if (not e->isLeaf())
                return false;

            /** \todo Why do I need those const_casts here? */
            if (refCount>=1)
                const_cast<FoamGridEntityImp<2,dimworld>*>(getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::REFINE;
            else if (refCount<0)
                const_cast<FoamGridEntityImp<2,dimworld>*>(getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::COARSEN;
            else
                const_cast<FoamGridEntityImp<2,dimworld>*>(getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::DO_NOTHING;

            return true;
        }
        
        /** \brief Return refinement mark for entity
        *
        * \return refinement mark (1,0,-1)
        */
        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            if (getRealImplementation(*e).target_->markState_ == FoamGridEntityImp<2,dimworld>::REFINE)
                return 1;
            if (getRealImplementation(*e).target_->markState_ == FoamGridEntityImp<2,dimworld>::COARSEN)
                return -1;
            
            return 0;
        }

        //! \brief Book-keeping routine to be called before adaptation
        bool preAdapt() {
            DUNE_THROW(NotImplemented, "preAdapt");
        }
        
        
        //! Triggers the grid refinement process
        bool adapt()
        {
            DUNE_THROW(NotImplemented, "adapt");
        }

        /** \brief Clean up refinement markers */
        void postAdapt() {
            DUNE_THROW(NotImplemented, "postAdapt");
        }
        
        /*@}*/
        
        /** @name Methods for parallel computations */
        /*@{*/
        
        /** \brief Size of the overlap on the leaf level */
        unsigned int overlapSize(int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the ghost cell layer on the leaf level */
        unsigned int ghostSize(int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the overlap on a given level */
        unsigned int overlapSize(int level, int codim) const {
            return 0;
        }
        
        
        /** \brief Size of the ghost cell layer on a given level */
        unsigned int ghostSize(int level, int codim) const {
            return 0;
        }
        
            
#if 0
        /** \brief Distributes this grid over the available nodes in a distributed machine
        *
        * \param minlevel The coarsest grid level that gets distributed
        * \param maxlevel does currently get ignored
        */
        void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
            DUNE_THROW(NotImplemented, "FoamGrid::loadBalance()");
        }
        
        /** \brief The communication interface
        *  @param T: array class holding data associated with the entities
        *  @param P: type used to gather/scatter data in and out of the message buffer
        *  @param codim: communicate entites of given codim
        *  @param if: one of the predifined interface types, throws error if it is not implemented
        *  @param level: communicate for entities on the given level
        *
        *  Implements a generic communication function sending an object of type P for each entity
        *  in the intersection of two processors. P has two methods gather and scatter that implement
        *  the protocol. Therefore P is called the "protocol class".
        */
        template<class T, template<class> class P, int codim>
        void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);
        
        /*! The new communication interface
        
        communicate objects for all codims on a given level
        */
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
        {}
        
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {}
#endif
        
#if 0        
        /** dummy collective communication */
        const CollectiveCommunication& comm () const
        {
            return ccobj_;
        }
#endif   
        /*@}*/
        
        
        // **********************************************************
        // End of Interface Methods
        // **********************************************************
        
    private:

        //! compute the grid indices and ids
    void setIndices()
    {
        // //////////////////////////////////////////
        //   Create the index sets
        // //////////////////////////////////////////
        for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
            FoamGridLevelIndexSet<const FoamGrid >* p
                = new FoamGridLevelIndexSet<const FoamGrid >();
            levelIndexSets_.push_back(p);
        }
        
        for (int i=0; i<=maxLevel(); i++)
            if (levelIndexSets_[i])
                levelIndexSets_[i]->update(*this, i);

        // Update the leaf indices
        leafIndexSet_.update(*this);
        
        // IdSets don't need updating

        }
#if 0         
        //! \todo Please doc me !
        CollectiveCommunication ccobj_;
#endif
    // Stores the lists of vertices, edges, elements for each level
    std::vector<tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                      std::list<FoamGridEntityImp<1,dimworld>>,
                      std::list<FoamGridEntityImp<2,dimworld>> > > entityImps_;

        //! Our set of level indices
        std::vector<FoamGridLevelIndexSet<const FoamGrid>*> levelIndexSets_;
        
        //! The leaf index set
        FoamGridLeafIndexSet<const FoamGrid > leafIndexSet_;
    
        //! The global id set
        FoamGridGlobalIdSet<const FoamGrid > globalIdSet_;
    
        //! The local id set
        FoamGridLocalIdSet<const FoamGrid > localIdSet_;

    /** \brief Counters that always provide the next free id for each dimension */
    array<unsigned int, dimension+1> freeIdCounter_;
    
}; // end Class FoamGrid




namespace Capabilities
{
    /** \brief True if the grid implements entities of a given codim.
      * 
      * FoamGrid implements all codimensions, hence this is always true
      */
    template<int dimworld,int codim>
    struct hasEntity< FoamGrid<dimworld>, codim>
    {
        static const bool v = true;
    };
    
    
    /** \brief True if the grid can be run on a distributed machine
      */
    template <int dimworld>
    struct isParallel< FoamGrid<dimworld> >
    {
        static const bool v = false;
    };
    
    
    //! \todo Please doc me !
    template<int dimworld>
    struct isLevelwiseConforming< FoamGrid<dimworld> >
    {
        static const bool v = true;
    };

    //! \todo Please doc me !
    template<int dimworld>
    struct isLeafwiseConforming< FoamGrid<dimworld> >
    {
        static const bool v = true;
    };
}

} // namespace Dune


// The factory should be present whenever the user includes foamgrid.hh.
// However since the factory needs to know the grid the include directive
// comes here at the end.
#include <dune/foamgrid/foamgrid/foamgridfactory.hh>

#endif
