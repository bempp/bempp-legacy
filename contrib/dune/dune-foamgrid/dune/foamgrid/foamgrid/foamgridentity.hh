#ifndef DUNE_FOAMGRID_ENTITY_HH
#define DUNE_FOAMGRID_ENTITY_HH

/** \file
* \brief The FoamGridEntity class
*/


namespace Dune {


// Forward declarations

template<int codim, int dim, class GridImp>
class FoamGridEntity;

template<int codim, class GridImp>
class FoamGridEntityPointer;

template<int codim, PartitionIteratorType pitype, class GridImp>
class FoamGridLevelIterator;

template<class GridImp>
class FoamGridLevelIntersectionIterator;

template<class GridImp>
class FoamGridLeafIntersectionIterator;

template<class GridImp>
class FoamGridHierarchicIterator;




template<int codim, int dim, class GridImp>
class FoamGridMakeableEntity :
    public GridImp::template Codim<codim>::Entity
{

    enum {dimworld = GridImp::dimensionworld};

    public:

        //! \todo Please doc me !
    FoamGridMakeableEntity(const FoamGridEntityImp<dim-codim,dimworld>* target) :
            GridImp::template Codim<codim>::Entity (FoamGridEntity<codim, dim, const GridImp>(target))
        {}
        
        
        //! \todo Please doc me !
    void setToTarget(const FoamGridEntityImp<dim-codim,dimworld>* target) {
            this->realEntity.setToTarget(target);
        }
        
        
        //! \todo Please doc me !
    const FoamGridEntityImp<dim-codim,dimworld>* getTarget() {
            return this->realEntity.target_;
        }
    
};


/** \brief The implementation of entities in a FoamGrid
*   \ingroup FoamGrid
*
*
*/
template<int codim, int dim, class GridImp>
class FoamGridEntity :
    public EntityDefaultImplementation <codim,dim,GridImp,FoamGridEntity>
{
    friend class FoamGridMakeableEntity<codim,dim,GridImp>;

    template <class GridImp_>
    friend class FoamGridLevelIndexSet;

    template <class GridImp_>
    friend class FoamGridLeafIndexSet;

    template <class GridImp_>
    friend class FoamGridLocalIdSet;

    template <class GridImp_>
    friend class FoamGridGlobalIdSet;

    friend class FoamGridEntityPointer<codim,GridImp>;

    
    private:
        
        typedef typename GridImp::ctype ctype;

    enum{dimworld = GridImp::dimensionworld};
        
    public:
    
        typedef typename GridImp::template Codim<codim>::Geometry Geometry;
    
        
        //! Constructor for an entity in a given grid level
    FoamGridEntity(const FoamGridEntityImp<dim-codim,GridImp::dimensionworld>* target) :
            target_(target)
        {}
        
        /** \brief Copy constructor */
        FoamGridEntity(const FoamGridEntity& original) :
            target_(original.target_)
        {}
    
    
        //! \todo Please doc me !
        FoamGridEntity& operator=(const FoamGridEntity& original)
        {
            if (this != &original)
            {
                geo_.reset();
                target_ = original.target_;
            }
            return *this;
        }


        //! level of this element
        int level () const {
            return target_->level();
        }
    
        
        /** \brief The partition type for parallel computing
        */
        PartitionType partitionType () const {
            return target_->partitionType();
        }
    
        
        /** Intra-element access to entities of codimension cc > codim. Return number of entities
        * with codimension cc.
        */
        template<int cc> int count () const{
            return target_->template count<cc>();
        }
        
        
        //! geometry of this entity
        const Geometry& geometry () const
        {
            if (geo_.get()==0)
                geo_ = std::auto_ptr<MakeableInterfaceObject<Geometry> >(new MakeableInterfaceObject<Geometry>(FoamGridGeometry<dim-codim,dimworld,GridImp>()));

            std::vector<FieldVector<double,dimworld> > coordinates(target_->corners());
            for (size_t i=0; i<target_->corners(); i++)
                coordinates[i] = target_->corner(i);

            GridImp::getRealImplementation(*geo_).setup(target_->type(), coordinates);
            return *geo_;
        }
    
    const FoamGridEntityImp<dim-codim,GridImp::dimensionworld>* target_;

        
    private:
    
        //! \todo Please doc me !
        void setToTarget(const FoamGridEntityImp<dim-codim,GridImp::dimensionworld>* target)
        {
            geo_.reset();
            target_ = target;
        }
    
        
        //! the current geometry
    mutable std::auto_ptr<MakeableInterfaceObject<Geometry> > geo_;
};




/** \brief Specialization for codim-0-entities, i.e., elements.
* \ingroup FoamGrid
*
* This class embodies the topological parts of elements of the grid.
* It has an extended interface compared to the general entity class.
* For example, Entities of codimension 0  allow to visit all neighbors.
*/
template<int dim, class GridImp>
class FoamGridEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, FoamGridEntity>
{

    enum {dimworld = GridImp::dimensionworld};

    public:
    
        typedef typename GridImp::template Codim<0>::Geometry Geometry;
    
        typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;
    
        //! The Iterator over intersections on this level
        typedef FoamGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;
    
        //! The Iterator over intersections on the leaf level
        typedef FoamGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;
    
        //! Iterator over descendants of the entity
        typedef FoamGridHierarchicIterator<GridImp> HierarchicIterator;
        
        
        //! Constructor for an entity in a given grid level
        FoamGridEntity(const FoamGridEntityImp<2,dimworld>* hostEntity) :
            geoInFather_(0),
            target_(hostEntity)
        {}
        
        
        /** \brief Copy constructor */
        FoamGridEntity(const FoamGridEntity& original) :
            geoInFather_(0),
            target_(original.target_)
        {}
    
        
        //! \todo Please doc me !
        FoamGridEntity& operator=(const FoamGridEntity& original)
        {
            if (this != &original)
            {
                geo_.reset();
                geoInFather_.reset();
                target_ = original.target_;
            }
            return *this;
        }
    
        
        //! Level of this element
        int level () const
        {
            return target_->level_;
        }
    
        
        /** \brief The partition type for parallel computing */
        PartitionType partitionType () const {
            return InteriorEntity;
        }
    
        
        //! Geometry of this entity
        const Geometry& geometry () const
        {
            if (not geo_.get())  // allocate on first use only
                geo_ = std::auto_ptr<MakeableInterfaceObject<Geometry> >(new MakeableInterfaceObject<Geometry>(FoamGridGeometry<dim,dimworld,GridImp>()));

            std::vector<FieldVector<double, dimworld> > coordinates(target_->corners());
            for (size_t i=0; i<target_->corners(); i++)
                coordinates[i] = target_->vertex_[i]->pos_;

            GridImp::getRealImplementation(*geo_).setup(target_->type(), coordinates);
            return *geo_;
        }
    
        
        /** \brief Return the number of subEntities of codimension cc.
        */
        template<int cc>
        int count () const
        {
            dune_static_assert(0<=cc && cc<=2, "Only codimensions with 0 <= cc <= 2 are valid!");
            return (cc==0) ? 1 : 3;
        }
        
    /** \brief Return level index of sub entity with codim = cc and local number i
     */
    int subLevelIndex (int i,unsigned int codim) const {
        assert(0<=codim && codim<=dim);
        switch (codim) {
        case 0:
            return target_->levelIndex_;
        case 1:
            return target_->edges_[i]->levelIndex_;
        case 2:
            return target_->vertex_[i]->levelIndex_;
        }
        DUNE_THROW(GridError, "Non-existing codimension requested!");
    }

    /** \brief Return leaf index of sub entity with codim = cc and local number i
     */
    int subLeafIndex (int i,unsigned int codim) const {
        assert(0<=codim && codim<=dim);
        switch (codim) {
        case 0:
            return target_->leafIndex_;
        case 1:
            return target_->edges_[i]->leafIndex_;
        case 2:
            return target_->vertex_[i]->leafIndex_;
        }
        DUNE_THROW(GridError, "Non-existing codimension requested!");
    }

    /** \brief Return index of sub entity with codim = cc and local number i
     */
    int subId (int i,unsigned int codim) const {
        assert(0<=codim && codim<=dim);
        switch (codim) {
        case 0:
            return target_->id_;
        case 1:
            return target_->edges_[i]->id_;
        case 2:
            return target_->vertex_[i]->id_;
        }
        DUNE_THROW(GridError, "Non-existing codimension requested!");
    }
    
        

        /** \brief Provide access to sub entity i of given codimension. Entities
        *  are numbered 0 ... count<cc>()-1
        */
        template<int codim>
        typename GridImp::template Codim<codim>::EntityPointer subEntity (int i) const{
            if (codim==0) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim,GridImp>( (FoamGridEntityImp<dim-codim,dimworld>*)this->target_);
            } else if (codim==1) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim,GridImp>( (FoamGridEntityImp<dim-codim,dimworld>*)this->target_->edges_[i]);
            } else if (codim==2) {
                // The cast is correct when this if clause is executed
                return FoamGridEntityPointer<codim,GridImp>( (FoamGridEntityImp<dim-codim,dimworld>*)this->target_->vertex_[i]);
            }
        }
    
        
        //! First level intersection
        FoamGridLevelIntersectionIterator<GridImp> ilevelbegin () const{
            return FoamGridLevelIntersectionIterator<GridImp>(target_, 0);
        }
    
        
        //! Reference to one past the last neighbor
        FoamGridLevelIntersectionIterator<GridImp> ilevelend () const{
            return FoamGridLevelIntersectionIterator<GridImp>(target_);
        }
    
        
        //! First leaf intersection
        FoamGridLeafIntersectionIterator<GridImp> ileafbegin () const{
            return FoamGridLeafIntersectionIterator<GridImp>(target_, (isLeaf()) ? 0 : target_->edges_.size());
        }
    
        
        //! Reference to one past the last leaf intersection
        FoamGridLeafIntersectionIterator<GridImp> ileafend () const{
            return FoamGridLeafIntersectionIterator<GridImp>(target_);
        }
    
        
        //! returns true if Entity has NO children
        bool isLeaf() const {
            return target_->isLeaf();
        }
    
    /** \brief Return true if this element has a father element */
    bool hasFather() const {
        return level()>0;
    }
        
        //! Inter-level access to father element on coarser grid.
        //! Assumes that meshes are nested.
        FoamGridEntityPointer<0,GridImp> father () const {
            return FoamGridEntityPointer<0,GridImp>(target_->father_);
        }
    
    
        /** \brief Location of this element relative to the reference element element of the father.
        * This is sufficient to interpolate all dofs in conforming case.
        * Nonconforming may require access to neighbors of father and
        * computations with local coordinates.
        * On the fly case is somewhat inefficient since dofs  are visited several times.
        * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
        * implementation of numerical algorithms is only done for simple discretizations.
        * Assumes that meshes are nested.
        */
        const LocalGeometry& geometryInFather () const {
            DUNE_THROW(NotImplemented, "geometryInFather");
//             if (geoInFather_==0)
//                 geoInFather_ = new MakeableInterfaceObject<LocalGeometry>(hostEntity_->geometryInFather());
            return *geoInFather_;
        }
    
        
        /** \brief Inter-level access to son elements on higher levels<=maxlevel.
        * This is provided for sparsely stored nested unstructured meshes.
        * Returns iterator to first son.
        */
        FoamGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
        {
            FoamGridHierarchicIterator<GridImp> it(maxLevel);

            // Load sons of old target onto the iterator stack
            if (level()<=maxLevel && !isLeaf())
                for (size_t i=0; i<target_->nSons(); i++)
                    it.elemStack.push(target_->sons_[i]);
            
            it.virtualEntity_.setToTarget((it.elemStack.empty()) 
                                          ? NULL : it.elemStack.top());
            
            return it;
        }
    
        
        //! Returns iterator to one past the last son
        FoamGridHierarchicIterator<GridImp> hend (int maxLevel) const
        {
            return FoamGridHierarchicIterator<const GridImp>(maxLevel);
        }
        
        
        // /////////////////////////////////////////
        //   Internal stuff
        // /////////////////////////////////////////
    
        
        /** \brief Make this class point to a new FoamGridEntityImp object */
        void setToTarget(const FoamGridEntityImp<2,dimworld>* target)
        {
            geo_.reset();
            geoInFather_.reset();
            target_ = target;
        }
        
        //! the current geometry
        mutable std::auto_ptr<MakeableInterfaceObject<Geometry> > geo_;
        
        /** \brief The geometry of this element as embedded in its father (if there is one) */
        mutable std::auto_ptr<MakeableInterfaceObject<LocalGeometry> > geoInFather_;

    const FoamGridEntityImp<2,dimworld>* target_;
        
    private:
    
        typedef typename GridImp::ctype ctype;

}; // end of FoamGridEntity codim = 0


} // namespace Dune


#endif
