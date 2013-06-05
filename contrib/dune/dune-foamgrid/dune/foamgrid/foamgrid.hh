// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_HH
#define DUNE_FOAMGRID_HH

/** \file
* \brief The FoamGrid class
*/

#include <list>
#include <set>

#include <dune/common/collectivecommunication.hh>
#include <dune/common/tuples.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// Implementation classes
#include "foamgrid/foamgridvertex.hh"
#include "foamgrid/foamgridedge.hh"
#include "foamgrid/foamgridelements.hh"

// The components of the FoamGrid interface
#include "foamgrid/foamgridgeometry.hh"
#include "foamgrid/foamgridlocalgeometry.hh"
#include "foamgrid/foamgridentity.hh"
#include "foamgrid/foamgridentitypointer.hh"
#include "foamgrid/foamgridentityseed.hh"
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
        CollectiveCommunication<Dune::FoamGrid<dimworld> > ,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        FoamGridEntitySeed /*,
                             FoamGridLocalGeometry*/
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
        : globalRefined(), numBoundarySegments_()
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
            //  DUNE_THROW(Dune::NotImplemented, "numBoundarySegments");
            return numBoundarySegments_;
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
        
        /** \brief Create EntityPointer from EnitySeed */
        template < class EntitySeed >
        typename Traits::template Codim<EntitySeed::codimension>::EntityPointer
            entityPointer(const EntitySeed& seed) const
        {
            typedef FoamGridEntityPointer<EntitySeed::codimension, const FoamGrid> EntityPointerImpl;
            typedef typename Traits::template Codim<EntitySeed::codimension>::EntityPointer EntityPointer;

            return EntityPointer(EntityPointerImpl(seed.getImplementationPointer()));
        }

        
        /** @name Grid Refinement Methods */
        /*@{*/
        
        
        /** global refinement
        */
        void globalRefine (int refCount)
        {
            willCoarsen=false;
            
            if(maxLevel()+refCount<0)
                DUNE_THROW(GridError, "Grid has only "<<maxLevel()<<" levels. Cannot do "
                           <<" globalRefine("<<refCount<<")");
            // The leafiterator is simply successively visiting all levels from fine to coarse
            // and just checking whether the isLeaf flag is set. Using it to identify
            // elements that need refinement, would produce and endless loop, as the newly 
            // added elements
            // would later be identified as elements that still need refinement.
            std::size_t oldLevels =entityImps_.size();
            typedef typename 
                std::vector<tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                                  std::list<FoamGridEntityImp<1,dimworld> >,
                                  std::list<FoamGridEntityImp<2,dimworld> > > >::reverse_iterator
                                  LevelIterator;
            // Allocate space for the new levels. Thus the rend iterator will not get
            // invalid due to newly added levels.
            entityImps_.reserve(oldLevels+refCount);
            levelIndexSets_.reserve(oldLevels+refCount);
            LevelIterator level = entityImps_.rbegin();
            // Add tuples for the new levels.
            for(int i=0; i < refCount; ++i){
                entityImps_.push_back(tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                                            std::list<FoamGridEntityImp<1,dimworld> >,
                                            std::list<FoamGridEntityImp<2,dimworld> > >());
                levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
            }
            if(refCount < 0){
                if(globalRefined+refCount<0)
                {
                    for(int i=refCount; i<0; ++i)
                    {
                        // Mark each leaf element for coarsening.
                        typedef typename Traits::template Codim<0>::LeafIterator Iterator;
                        for(Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
                            elem != end; ++elem)
                        {
                            mark(-1,*elem);
                        }
                        // do the adaptation
                        preAdapt();
                        adapt();
                        postAdapt();
                    }
                    globalRefined=0;
                    return;
                }
                else
                {
                    
                
                    for(int i=refCount; i<0; ++i){
                        delete levelIndexSets_.back();
                        levelIndexSets_.pop_back();
                    }
                    entityImps_.resize(oldLevels+refCount);
                    
                    //To be able to create the leaf level we need to set
                    // the sons of the entities of maxlevel to null
                    typename std::list<FoamGridEntityImp<0,dimworld> >::iterator vIt    = 
                        Dune::get<0>(entityImps_[maxLevel()]).begin();
                    typename std::list<FoamGridEntityImp<0,dimworld> >::iterator vEndIt = 
                        Dune::get<0>(entityImps_[maxLevel()]).end();
                    for (; vIt!=vEndIt; ++vIt) {
                        vIt->son_=nullptr;
                    }
                    typename std::list<FoamGridEntityImp<1,dimworld> >::iterator edIt    = 
                        Dune::get<1>(entityImps_[maxLevel()]).begin();
                    typename std::list<FoamGridEntityImp<1,dimworld> >::iterator edEndIt = 
                        Dune::get<1>(entityImps_[maxLevel()]).end();
                    for (; edIt!=edEndIt; ++edIt) {
                        edIt->sons_[0]=nullptr;
                        edIt->sons_[1]=nullptr;
                        edIt->nSons_=0;
                    }
                    typename std::list<FoamGridEntityImp<2,dimworld> >::iterator elIt    = 
                        Dune::get<2>(entityImps_[maxLevel()]).begin();
                    typename std::list<FoamGridEntityImp<2,dimworld> >::iterator elEndIt = 
                        Dune::get<2>(entityImps_[maxLevel()]).end();
                    for (; elIt!=elEndIt; ++elIt) {
                        elIt->sons_[0]=nullptr;
                        elIt->sons_[1]=nullptr;
                        elIt->sons_[2]=nullptr;
                        elIt->sons_[3]=nullptr;
                        elIt->nSons_=0;
                    }
                }
                
            }else{
                
                // sanity check whether the above asssumption really holds.
                assert(&entityImps_[oldLevels-1]==&(*level));
                
                // Do the actual refinement.
                // We start with the finest level
                std::size_t levelIndex;
                for(levelIndex=oldLevels-1; level!=entityImps_.rend(); 
                    ++level, --levelIndex){
                    
                    typedef typename std::list<FoamGridEntityImp<2,dimworld> >::iterator ElementIterator;
                    bool foundLeaf=false;
                    
                    for(ElementIterator element=Dune::get<2>(*level).begin(); element != Dune::get<2>(*level).end(); ++element)
                        if(element->isLeaf()){
                            foundLeaf = true;
                            dverb<<"refining element "<< &(*element)<<std::endl;
                            if(element->type().isTriangle())
                                refineSimplexElement(*element, refCount);
                            else
                                DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
                        }
                    if(!foundLeaf)
                        break;
                    }
                for(levelIndex+=2;levelIndex!=entityImps_.size(); ++levelIndex)
                    levelIndexSets_[levelIndex]->update(*this, levelIndex);
            }
            
            // Update the leaf indices
            leafIndexSet_.update(*this);

            globalRefined=std::max(globalRefined+refCount,0);
            postAdapt();
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
                const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::REFINE;
            else if (refCount<0)
                const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::COARSEN;
            else
                const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*e).target_)->markState_ = FoamGridEntityImp<2,dimworld>::DO_NOTHING;

            return true;
        }
        
        /** \brief Return refinement mark for entity
        *
        * \return refinement mark (1,0,-1)
        */
        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            if (this->getRealImplementation(*e).target_->markState_ == FoamGridEntityImp<2,dimworld>::REFINE)
                return 1;
            if (this->getRealImplementation(*e).target_->markState_ == FoamGridEntityImp<2,dimworld>::COARSEN)
                return -1;
            
            return 0;
        }

        //! \brief Book-keeping routine to be called before adaptation
        bool preAdapt() {
            // Loop over all leaf entities and check whether they might be
            // coarsened. If there is one return true.
            typedef typename Traits::template Codim<0>::LeafIterator Iterator;
            int addLevels = 0;
            willCoarsen = false;
            
            for(Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
                elem != end; ++elem)
            {
                int mark=getMark(*elem);
                addLevels=std::max(addLevels, elem->level()+mark-maxLevel());
                
                if(mark<0)
                {
                    // If this element is marked for coarsening, but another child
                    // of this element's father is marked for refinement, then we 
                    // need to reset the marker to doNothing
                    bool otherChildRefined=false;
                    FoamGridEntityImp<2,dimworld>& father = *this->getRealImplementation(*elem).target_->father_;
                    typedef typename array<FoamGridEntityImp<2,dimworld>*,4>::iterator
                        ChildrenIter;
                    for(ChildrenIter child=father.sons_.begin(); child != 
                            father.sons_.end(); ++child)
                        otherChildRefined=otherChildRefined || 
                            (*child)->markState_==FoamGridEntityImp<2,dimworld>::REFINE;
                    
                    if(otherChildRefined)
                    {
                        for(ChildrenIter child=father.sons_.begin(); child != 
                                father.sons_.end(); ++child)
                            if((*child)->markState_==FoamGridEntityImp<2,dimworld>::COARSEN)
                                (*child)->markState_=FoamGridEntityImp<2,dimworld>::DO_NOTHING;
                    }
                    else
                        willCoarsen = willCoarsen || mark<0;
                }
            }
            if(addLevels)
            {
                entityImps_.push_back(tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                                            std::list<FoamGridEntityImp<1,dimworld> >,
                                            std::list<FoamGridEntityImp<2,dimworld> > >());
                levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
            }
            return willCoarsen;
        }
        
        
        //! Triggers the grid refinement process
        bool adapt()
        {
            std::set<std::size_t> levelsChanged;
            bool haveRefined=false;
            
            // Loop over all leaf elements and refine/coarsen those that marked for it.
            typedef typename Traits::template Codim<0>::LeafIterator Iterator;
            for(Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
                elem != end; ++elem)
            {
                int mark=getMark(*elem);
                if(mark>0)
                {
                    // Refine Triangles
                    if(elem->type().isTriangle())
                    {
                        refineSimplexElement(*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_), 1);
                        levelsChanged.insert(elem->level()+1);
                        haveRefined=true;
                    }
                    else
                        DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
                }
                if(mark<0) // If simplex was allready treated by coarsenSimplex mark will be 0
                {
                    // Coarsen triangles
                    if(elem->type().isTriangle())
                    {
                        assert(elem->level());
                        coarsenSimplexElement(*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_));
                                              levelsChanged.insert(elem->level());
                    }
                    else
                        DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
                }
            }
            if(!willCoarsen)
                return haveRefined;

            typedef typename std::set<std::size_t>::const_reverse_iterator SIter;
            for(SIter level=levelsChanged.rbegin(); level!=levelsChanged.rend(); ++level)
            {
                // First delete the pointer
                erasePointersToEntities(Dune::get<2>(entityImps_[*level]));
                // Now delete the actual vertices.
                eraseVanishedEntities(Dune::get<0>(entityImps_[*level]));
                // erase vanished edges
                {
                    typedef typename std::list<FoamGridEntityImp<1,dimworld> >::iterator
                        EdgeIter;
                    for(EdgeIter edge=Dune::get<1>(entityImps_[*level]).begin(),
                            edgeEnd=Dune::get<1>(entityImps_[*level]).end();
                        edge != Dune::get<1>(entityImps_[*level]).end();)
                    {
                        bool hasSameLevelElements=false;
                        typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator
                            ElementIter;
                        for(ElementIter elem=edge->elements_.begin();
                            elem != edge->elements_.end(); ++elem)
                            hasSameLevelElements = hasSameLevelElements || (*elem)->level()==*level;
                        assert(edge->willVanish_!=hasSameLevelElements);
                        
                        if(!hasSameLevelElements)
                        {
                            assert(edge->willVanish_);
                            // erase returns next position
                            edge=Dune::get<1>(entityImps_[*level]).erase(edge);
                        }else
                            // increment
                            ++edge;
                    }
                }
                // And the elements
                eraseVanishedEntities(Dune::get<2>(entityImps_[*level]));
                
                if( Dune::get<0>(entityImps_[*level]).size())
                {
                    assert(Dune::get<1>(entityImps_[*level]).size() &&
                           Dune::get<2>(entityImps_[*level]).size());
                    // Update the level indices.
                    levelIndexSets_[*level]->update(*this, *level);
                }
                else
                {
                    assert(!Dune::get<1>(entityImps_[*level]).size() &&
                           !Dune::get<2>(entityImps_[*level]).size());
                    if(static_cast<int>(*level)==maxLevel())
                        entityImps_.pop_back();
                }
            }
                                                           
            if(levelsChanged.size())
                // Update the leaf indices
            leafIndexSet_.update(*this);
            globalRefined=0;
            return haveRefined;
        }

        

        /** \brief Clean up refinement markers */
        void postAdapt() {
            willCoarsen=false;
            // Loop over all leaf entities and remove the isNew Marker.
            typedef typename Traits::template Codim<0>::LeafIterator Iterator;
            
            for(Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
                elem != end; ++elem)
            {
                FoamGridEntityImp<2,dimworld>& element=*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_);
                element.isNew_=false;
                assert(!element.willVanish_);
                if(element.father_)
                    element.father_->markState_=FoamGridEntityImp<2,dimworld>::DO_NOTHING;
            }
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
        
        /** dummy collective communication */
        const typename Traits::CollectiveCommunication& comm () const
        {
            return ccobj_;
        }
        /*@}*/
        
        
        // **********************************************************
        // End of Interface Methods
        // **********************************************************
        
    private:

        //! \brief erases pointers in father elements to vanished entities of the element
        void erasePointersToEntities(std::list<FoamGridEntityImp<2,dimworld> >& elements)
        {
            typedef typename std::list<FoamGridEntityImp<2,dimworld> >::iterator EntityIterator;
            for(EntityIterator element=elements.begin(); 
                element != elements.end(); ++element)
            {
                if(element->willVanish_)
                {
                    FoamGridEntityImp<2,dimworld>& father=*element->father_;
                    /*if(father.sons_[0]==nullptr)
                    {
                        assert(father.isLeaf());
                        // father was already processes for another element
                        return;
                        }*/
                    
                    for(unsigned int i=0; i<father.nSons_; i++)
                        father.sons_[i]=nullptr;
                    for(unsigned int i=0; i<father.corners(); i++)
                        if(father.vertex_[i]->son_!=nullptr)
                            if(father.vertex_[i]->son_->willVanish_)
                                father.vertex_[i]->son_=nullptr;
                    for(unsigned int i=0; i<father.corners(); i++)
                        if(father.edges_[i]->sons_[0]!=nullptr
                           && father.edges_[i]->sons_[0]->willVanish_)
                            for(unsigned int j=0; j<2; j++)
                            {
                                assert(father.edges_[i]->sons_[j]!=nullptr);
                                assert(father.edges_[i]->sons_[j]->willVanish_);
                                father.edges_[i]->sons_[j]=nullptr;
                            }        
                }
            }
            
        }
    
        //! \brief Erase Entities from memory that vanished due to coarsening.
        //! \warning This method has to be called first for i=0.
        //! \tparam i The dimension of the entities.
        //! \param  levelEntities The vector with the level entitied
        template<int i>
        void eraseVanishedEntities(std::list<FoamGridEntityImp<i,dimworld> >& levelEntities)
        {
            typedef typename std::list<FoamGridEntityImp<i,dimworld> >::iterator EntityIterator;
            for(EntityIterator entity=levelEntities.begin(); 
                entity != levelEntities.end();)
            {
                if(entity->willVanish_)
                {
                    entity=levelEntities.erase(entity);
                 }else
                    ++entity;
            }
        }


    //! \brief Coarsen an Element
    //! \param element The element to coarsen
    void coarsenSimplexElement(FoamGridEntityImp<2,dimworld>& element)
    {
        // If we coarsen an element, this means that we erase all chidren of its father
        // to prevent inconsistencies.
        const FoamGridEntityImp<2,dimworld>& father = *(element.father_);

        // The edges that might be erased
        std::set<FoamGridEntityImp<1,dimworld>*> childEdges;
        // The vertices that might be erased
        std::set<FoamGridEntityImp<0,dimworld>*> childVertices;
        typedef typename array<FoamGridEntityImp<2,dimworld>*,4>::const_iterator
            ChildrenIter;
        for(ChildrenIter child=father.sons_.begin(); child != 
                father.sons_.end(); ++child)
        {
            // Remember element for the actual deletion taking place later
            (*child)->markState_=FoamGridEntityImp<2,dimworld>::IS_COARSENED;
            (*child)->willVanish_=true;

            typedef typename array<FoamGridEntityImp<1,dimworld>*,3>::iterator
                EdgeIter;
            for(EdgeIter edge=(*child)->edges_.begin(); edge != (*child)->edges_.end();
                ++edge)
            {
                // Remove references to elements that will be erased
                typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
                    ElementIter;
                for(ElementIter element = (*edge)->elements_.begin(); 
                    element != (*edge)->elements_.end();
                    ++element)
                {
                    for(ChildrenIter child1=father.sons_.begin(); child1 != 
                            father.sons_.end(); ++child1)
                        if(*element==*child1)
                        {
                            // To prevent removing an element from a vector
                            // we just overwrite the element with its father.
                            // later we will check the edges and erase all
                            // of them that do not have any elements on the same
                            // level.
                            *element=&father;
                        }
                }
            }
                

                
            typedef FoamGridEntityImp<0,dimworld>** VertexIter;
            for(VertexIter vertex=(*child)->vertex_; vertex != (*child)->vertex_+(*child)->corners();
                ++vertex)
                childVertices.insert(*vertex);
            typedef typename array<FoamGridEntityImp<1,dimworld>*, 3>::iterator EdgeIter;
            for(EdgeIter edge=(*child)->edges_.begin(); edge!= (*child)->edges_.end();
                ++edge)
                childEdges.insert(*edge);
        }
        
        // Check whether those guys are really erased.
        // That is, we remove all vertices and edges that are part of a child element of one
        // of the neighbours of the father element
        typedef  typename array<FoamGridEntityImp<1,dimworld>*,3>::const_iterator
                EdgeIter;
        for(EdgeIter edge=father.edges_.begin(); edge != father.edges_.end();
            ++edge)
        {
            typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
                NeighborIter;
            for(NeighborIter neighbor = (*edge)->elements_.begin(); 
                neighbor != (*edge)->elements_.end();
                ++neighbor)
            {
                assert((*neighbor)->level()<=father.level());
                if(*neighbor != &father && (*neighbor)->level()==father.level() &&
                   !(*neighbor)->isLeaf())
                {
                    // This is a real neighbor element on the same level
                    // Check whether one of its children is marked for coarsening
                    bool coarsened=false;
                    for(ChildrenIter child=(*neighbor)->sons_.begin();
                        child != (*neighbor)->sons_.end(); ++child)
                        if((*child)->mightVanish())
                        {
                            coarsened=true;
                            break;
                        }
                        
                    if(!coarsened)
                    {
                        // Remove all entities that exist in the children of the neighbor
                        // from the list
                        for(ChildrenIter child=(*neighbor)->sons_.begin();
                            child != (*neighbor)->sons_.end(); ++child)
                        {
                            typedef typename array<FoamGridEntityImp<1,dimworld>*,3>::iterator
                                EdgeIter;
                            for(EdgeIter edge=(*child)->edges_.begin(); edge != (*child)->edges_.end();
                                ++edge)
                                childEdges.erase(*edge);
                            typedef FoamGridEntityImp<0,dimworld>** VertexIter;
                            for(VertexIter vertex=(*child)->vertex_; 
                                vertex != (*child)->vertex_+(*child)->corners();
                                ++vertex)
                                childVertices.erase(*vertex);
                        }
                    }
                }
            }
        }
        typedef typename std::set<FoamGridEntityImp<1,dimworld>*>::iterator SEdgeIter;
        for(SEdgeIter e=childEdges.begin(); e!=childEdges.end(); ++e)
            (*e)->willVanish_=true;
        typedef typename std::set<FoamGridEntityImp<0,dimworld>*>::iterator VertexIter;
        for(VertexIter v=childVertices.begin(); v!=childVertices.end(); ++v)
            (*v)->willVanish_=true;
    }
    
                        
                             
        

    //! \brief refine an Element
    //! \param elemet The element to refine
    //! \param refCount How many times to refine the element
    void refineSimplexElement(FoamGridEntityImp<2,dimworld>& element,
                       int refCount)
    {
        if(refCount<0)
        {
            // We always remove all the children from the father.
            // Removing means:
            // 1. remove pointers in father element, edges and vertices
            // 2. Mark removed entities for deletion.
            DUNE_THROW(NotImplemented, "Coarsening not implemented yet");
            return;
        }
        

        // TODO: Currently we are assuming that only globalRefine is available
        // and therefore the next level is not yet present.
        // For real adaptivity some of the edges and vertices might already present.
        // Therefore we need some kind of detection for this later on.

        unsigned int nextLevel=element.level()+1;

        
        std::cout<<"Vertices "<<nextLevel<<": "<<Dune::get<0>(entityImps_[nextLevel]).size()
                 <<std::endl;
        std::cout<<"Edges "<<nextLevel<<": "<<Dune::get<1>(entityImps_[nextLevel]).size()
                 <<std::endl;
        std::cout<<"Elements "<<nextLevel<<": "<<Dune::get<2>(entityImps_[nextLevel]).size()
                 <<std::endl;

        array<FoamGridEntityImp<0,dimworld>*, 6> nextLevelVertices;
        std::size_t vertexIndex=0;
        // create copies of the vertices of the element
        for(unsigned int c=0; c<element.corners(); ++c){
            dverb<<"Processing vertex "<<element.vertex_[c]<<std::endl;
            if(element.vertex_[c]->son_==nullptr){
                // Not refined yet
                Dune::get<0>(entityImps_[nextLevel])
                    .push_back(FoamGridEntityImp<0,dimworld>(nextLevel, 
                                                             element.vertex_[c]->pos_,
                                                             element.vertex_[c]->id_));
                FoamGridEntityImp<0,dimworld>& newVertex = 
                               Dune::get<0>(entityImps_[nextLevel]).back();
                element.vertex_[c]->son_=&newVertex;                
            }
            check_for_duplicates(nextLevelVertices, element.vertex_[c]->son_, vertexIndex);
            nextLevelVertices[vertexIndex++]=element.vertex_[c]->son_;
        }
        
        // create new vertices from edge-midpoints together with the new edges that 
        // have a father
        typedef typename array<FoamGridEntityImp<1,dimworld>*, 3>::iterator EdgeIterator;
        
        array<FoamGridEntityImp<1,dimworld>*, 9> nextLevelEdges;
        std::size_t edgeIndex=0;
        const Dune::GenericReferenceElement<double,dimension>& refElement
            = Dune::GenericReferenceElements<double, dimension>::general(element.type());

        // I am just to dumb for a general edge to verticex mapping.
        // Therefore we just store it here
        array<std::pair<unsigned int,unsigned int>,3 > edgeVertexMapping;
        edgeVertexMapping[0]=std::make_pair(0,1);
        edgeVertexMapping[1]=std::make_pair(2,0);
        edgeVertexMapping[2]=std::make_pair(1,2);
        
        for(EdgeIterator edge=element.edges_.begin(); edge != element.edges_.end(); ++edge){
            typedef FoamGridEntityImp<0,dimworld> FoamGridVertex;
            const FoamGridVertex* v0 = element.vertex_[refElement.subEntity(edgeIndex/2, 1, 0, 2)];
            const FoamGridVertex* v1 = element.vertex_[refElement.subEntity(edgeIndex/2, 1, 1, 2)];
            if(!(*edge)->nSons_){
                // Not refined yet
                // Compute edge midpoint
                FieldVector<double, dimworld> midPoint;
                for(int dim=0; dim<dimworld;++dim)
                    midPoint[dim]=((*edge)->vertex_[0]->pos_[dim] 
                                   + (*edge)->vertex_[1]->pos_[dim]) /2.0;
                
                //create midpoint
                Dune::get<0>(entityImps_[nextLevel])
                    .push_back(FoamGridEntityImp<0,dimworld>(nextLevel, midPoint, 
                                                             freeIdCounter_[0]++));
                FoamGridEntityImp<0,dimworld>& midVertex = 
                    Dune::get<0>(entityImps_[nextLevel]).back();
                check_for_duplicates(nextLevelVertices, &midVertex, vertexIndex);
                nextLevelVertices[vertexIndex++]=&midVertex;

                // sanity check for DUNE numbering
                dvverb<<"edge "<<edgeIndex/2<<": "<<"("<<(*edge)->vertex_[0]->son_<<","<<(*edge)->vertex_[1]->son_<<") with father ("<<(*edge)->vertex_[0]<<","<<(*edge)->vertex_[1]<<")"<<std::endl;
                assert(v0->son_!=nullptr);
                assert(v1->son_!=nullptr);
                assert(v0->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].first] ||
                       v0->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].second]);

                assert(v1->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].first] ||
                       v1->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].second]);

                // create the edges and publish them in the father
                Dune::get<1>(entityImps_[nextLevel])
                    .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[edgeVertexMapping[edgeIndex/2].first], &midVertex, 
                                                             nextLevel, freeIdCounter_[1]++, *edge));
                (*edge)->sons_[0] = &Dune::get<1>(entityImps_[nextLevel]).back();
                ++((*edge)->nSons_);
                // Inherit the boundaryId
                (*edge)->sons_[0]->boundaryId_=(*edge)->boundaryId_;
                nextLevelEdges[edgeIndex++]= (*edge)->sons_[0];

                // Initialize the elements_ vector of the new edge
                // with that of the father. Later we will overwrite it
                // with the correct values.
                (*edge)->sons_[0]->elements_=(*edge)->elements_;

                assert((*edge)->vertex_[1]->son_!=nullptr);
                Dune::get<1>(entityImps_[nextLevel])
                    .push_back(FoamGridEntityImp<1,dimworld>(&midVertex, nextLevelVertices[edgeVertexMapping[edgeIndex/2].second], 
                                                             nextLevel, freeIdCounter_[1]++, *edge));
                (*edge)->sons_[1] = &Dune::get<1>(entityImps_[nextLevel]).back();
                ++((*edge)->nSons_);
                // Inherit the boundaryId
                (*edge)->sons_[1]->boundaryId_=(*edge)->boundaryId_;
                nextLevelEdges[edgeIndex++]= (*edge)->sons_[1];
                
                // Initialize the elements_ vector of the new edge
                // with that of the father. Later we will overwrite it
                // with the correct values.
                (*edge)->sons_[1]->elements_=(*edge)->elements_;
        
                (*edge)->nSons_=2;
            }else{
                // Edges do already exist. Just add its sons to nextLevelEdges
                // but make sure that the one containing vertex edgeIndex comes first
                if((*edge)->sons_[0]->vertex_[0]->id_ == 
                   nextLevelVertices[edgeVertexMapping[edgeIndex/2].first]->id_ ||
                   (*edge)->sons_[0]->vertex_[1]->id_ == 
                   nextLevelVertices[edgeVertexMapping[edgeIndex/2].first]->id_){
                    nextLevelEdges[edgeIndex++]=(*edge)->sons_[0];
                    nextLevelEdges[edgeIndex++]=(*edge)->sons_[1];
                }else{
                    nextLevelEdges[edgeIndex++]=(*edge)->sons_[1];
                    nextLevelEdges[edgeIndex++]=(*edge)->sons_[0];
                }
                if((*edge)->sons_[0]->vertex_[0]->id_!=(*edge)->vertex_[0]->id_ &&
                   (*edge)->sons_[0]->vertex_[0]->id_!=(*edge)->vertex_[1]->id_){
                    //vertex 0 is the midpoint
                    check_for_duplicates(nextLevelVertices, (*edge)->sons_[0]->vertex_[0], vertexIndex);
                    nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0,dimworld>*>((*edge)->sons_[0]->vertex_[0]);
                }else{
                    check_for_duplicates(nextLevelVertices, (*edge)->sons_[0]->vertex_[1], vertexIndex);
                    nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0,dimworld>*>((*edge)->sons_[0]->vertex_[1]);
                }
            }
        }
        assert(edgeIndex==6);
        // Create the edges that lie within the father element
        // first the one that lies opposite to the vertex 0 in the father
        Dune::get<1>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[3],
                                                     nextLevelVertices[4], nextLevel, 
                                                     freeIdCounter_[1]++));
        nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();
        
        // the one opposite to father vertex 1
        Dune::get<1>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[3],
                                                     nextLevelVertices[5], nextLevel, 
                                                     freeIdCounter_[1]++));
        nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();
        
        // and the one opposite to father vertex 2
        Dune::get<1>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[4],
                                                     nextLevelVertices[5], nextLevel, 
                                                     freeIdCounter_[1]++));
        nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();
        
        assert(edgeIndex==nextLevelEdges.size());
        assert(vertexIndex==nextLevelVertices.size());
        
        array<FoamGridEntityImp<2,dimworld>*, 4> nextLevelElements;
        // create the new triangles that lie in the corners
        // First the one that contains vertex 0 of the father.
        Dune::get<2>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
            
        FoamGridEntityImp<2,dimworld>* newElement = &(Dune::get<2>(entityImps_[nextLevel])
                                                      .back());
        newElement->isNew_=true;
        newElement->father_=&element;
        newElement->edges_[0]=nextLevelEdges[0];
        newElement->edges_[1]=nextLevelEdges[3];
        newElement->edges_[2]=nextLevelEdges[6];
        newElement->vertex_[0]=nextLevelVertices[0];
        newElement->vertex_[1]=nextLevelVertices[3];
        newElement->vertex_[2]=nextLevelVertices[4];
        newElement->refinementIndex_=0;
        nextLevelElements[0]=newElement;
        element.sons_[0]=newElement;
        dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;
        
        // Next the one that contains vertex 1 of the father.
        Dune::get<2>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
        newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
        newElement->isNew_=true;
        newElement->father_=&element;
        newElement->edges_[0]=nextLevelEdges[4];
        newElement->edges_[1]=nextLevelEdges[1];
        newElement->edges_[2]=nextLevelEdges[7];
        newElement->vertex_[0]=nextLevelVertices[1];
        newElement->vertex_[1]=nextLevelVertices[5];
        newElement->vertex_[2]=nextLevelVertices[3];
        newElement->refinementIndex_=1;
        nextLevelElements[1]=newElement;
        element.sons_[1]=newElement;
        dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


        // Last the one that contains vertex 2 of the father.
        Dune::get<2>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
        newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
        newElement->isNew_=true;
        newElement->father_=&element;
        newElement->edges_[0]=nextLevelEdges[2];
        newElement->edges_[1]=nextLevelEdges[5];
        newElement->edges_[2]=nextLevelEdges[8];
        newElement->vertex_[0]=nextLevelVertices[2];
        newElement->vertex_[1]=nextLevelVertices[4];
        newElement->vertex_[2]=nextLevelVertices[5];
        newElement->refinementIndex_=2;
        nextLevelElements[2]=newElement;
        element.sons_[2]=newElement;
        dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


        // create the triangle in the center
        Dune::get<2>(entityImps_[nextLevel])
            .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
        newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
        newElement->isNew_=true;
        newElement->father_=&element;
        newElement->edges_[0]=nextLevelEdges[7];
        newElement->edges_[1]=nextLevelEdges[6];
        newElement->edges_[2]=nextLevelEdges[8];
        newElement->vertex_[0]=nextLevelVertices[3];
        newElement->vertex_[1]=nextLevelVertices[5];
        newElement->vertex_[2]=nextLevelVertices[4];
        newElement->refinementIndex_=3;
        nextLevelElements[3]=newElement;
        element.sons_[3]=newElement;
        dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


        // Now that all the triangle are created, we can update the elements attached
        // to the edges.
        // The new neighbors of the edges lying on edges of the father element.
        std::size_t neighbors[6] = {0, 1, 2, 0, 1, 2};
        dvverb<<" element "<<&element<<std::endl;
        for(std::size_t i=0; i<6; ++i){            
            // Overwrite the father element by the newly created elements.
            dvverb<<" neighbour "<<i<<": ";
            overwriteFineLevelNeighbours(*nextLevelEdges[i], nextLevelElements[neighbors[i]],
                                         &element);
        }

        // Update the neighbours of the inner edges
        nextLevelEdges[6]->elements_.push_back(nextLevelElements[0]);
        nextLevelEdges[6]->elements_.push_back(nextLevelElements[3]);
        nextLevelEdges[7]->elements_.push_back(nextLevelElements[3]);
        nextLevelEdges[7]->elements_.push_back(nextLevelElements[1]);
        nextLevelEdges[8]->elements_.push_back(nextLevelElements[3]);
        nextLevelEdges[8]->elements_.push_back(nextLevelElements[2]);
#ifndef NDEBUG
        for(std::size_t vidx=0; vidx<4; ++vidx){
            dvverb<<std::endl<<" Element "<<nextLevelElements[vidx]<<": ";
            for(EdgeIterator edge=nextLevelElements[vidx]->edges_.begin(); 
                edge != nextLevelElements[vidx]->edges_.end(); ++edge){
                dvverb<<std::endl<<"   edge "<<*edge<<": ";
                typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
                    ElementIterator;
                bool selfFound=false;
                for(ElementIterator elem=(*edge)->elements_.begin(); 
                    elem != (*edge)->elements_.end();
                    ++elem){
                    dvverb << *elem<<" ";
                    if(*elem == nextLevelElements[vidx])
                        selfFound=true;
                }
                assert(selfFound);
            }
        }
        
#endif
        element.nSons_=4;
        
        if((refCount--)>1){
            int i=0;
            typedef typename array<FoamGridEntityImp<2,dimworld>*, 4>::iterator
                ElementIterator;
            for(ElementIterator elem=nextLevelElements.begin();
                elem != nextLevelElements.end(); ++elem){
                dvverb<<std::endl<<"Refinining "<<(*elem)<<" (son of"<<&element<<") refCount="<<refCount<<" child="<<i++<<std::endl;
                refineSimplexElement(**elem, refCount);
            }
        }
        std::cout<<"end refineSimplex"<<std::endl;
        std::cout<<"Vertices "<<nextLevel<<": "<<Dune::get<0>(entityImps_[nextLevel]).size()
                 <<std::endl;
        std::cout<<"Edges "<<nextLevel<<": "<<Dune::get<1>(entityImps_[nextLevel]).size()
                 <<std::endl;
        std::cout<<"Elements "<<nextLevel<<": "<<Dune::get<2>(entityImps_[nextLevel]).size()
                 <<std::endl;
    }

    /**
     * \brief Overwrites the nighbours of this and descendant edges
     *
     * After returning all neighbours the previously pointed to the 
     * father will point to the son element
     * \param edge The edge to start overwriting with.
     * \param son The son element to substitute the father with.
     * \param father Pointer to the father element that is to be subsitituted.
     */
    void overwriteFineLevelNeighbours(FoamGridEntityImp<1,dimworld>& edge,
                                      FoamGridEntityImp<2,dimworld>* son,
                                      FoamGridEntityImp<2,dimworld>* father){
        typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
                ElementIterator;
#ifndef NDEBUG
        bool fatherFound=false;
#endif            
        for(ElementIterator elem=edge.elements_.begin(); 
            elem != edge.elements_.end();
            ++elem){
            dvverb << *elem<<" ";
            if(*elem == father){
#ifndef NDEBUG
                fatherFound=true;
#endif
                *elem = son;
            }
        }

#ifndef NDEBUG
        dvverb<<std::endl;
        assert(fatherFound);
#endif
        for(std::size_t i=0; i<edge.nSons_; ++i)
            overwriteFineLevelNeighbours(*edge.sons_[i], son, father);    
    }
    
    template<class C, class T>
    void check_for_duplicates(C& cont, const T& elem, std::size_t vertexIndex)
    {
#ifndef NDEBUG
        for(std::size_t i=0; i<vertexIndex; ++i)
            assert(cont[i]!=elem);
#endif
    }
    
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

        //! Collective communication interface
        typename Traits::CollectiveCommunication ccobj_;

    // Stores the lists of vertices, edges, elements for each level
    std::vector<tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                      std::list<FoamGridEntityImp<1,dimworld> >,
                      std::list<FoamGridEntityImp<2,dimworld> > > > entityImps_;

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

    /** \brief How many times was the leaf level globally refined. */
    int globalRefined;

    /** \brief The number of boundary segements. */
    std::size_t numBoundarySegments_;

    // True if the last call to preadapt returned true
    bool willCoarsen;
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
        static const bool v = false;
    };

    //! \todo Please doc me !
    template<int dimworld>
    struct isLeafwiseConforming< FoamGrid<dimworld> >
    {
        static const bool v = false;
    };
}

} // namespace Dune


// The factory should be present whenever the user includes foamgrid.hh.
// However since the factory needs to know the grid the include directive
// comes here at the end.
#include <dune/foamgrid/foamgrid/foamgridfactory.hh>

#endif
