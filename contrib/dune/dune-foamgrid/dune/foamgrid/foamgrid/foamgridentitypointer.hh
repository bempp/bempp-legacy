#ifndef DUNE_FOAMGRID_ENTITY_POINTER_HH
#define DUNE_FOAMGRID_ENTITY_POINTER_HH

/** \file
* \brief The FoamGridEntityPointer class
*/

namespace Dune {


/** Acts as a pointer to an  entities of a given codimension.
*/
template<int codim, class GridImp>
class FoamGridEntityPointer
{
    private:
    
    enum { dim      = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };
    
    public:
    
    //! export the type of the EntityPointer Implementation.
    //! Necessary for the typeconversion between Iterators and EntityPointer
    typedef FoamGridEntityPointer EntityPointerImp;

    /** \brief Codimension of entity pointed to */
    enum { codimension = codim };

        typedef typename GridImp::template Codim<codim>::Entity Entity;
        
    //! Constructor from a FoamGrid entity
    FoamGridEntityPointer (const FoamGridEntity<codim,dim,GridImp>& entity)
        : virtualEntity_(entity.target_) 
    {}

    FoamGridEntityPointer (const typename std::list<FoamGridEntityImp<dim-codim,dimworld> >::const_iterator& it)
        : virtualEntity_(&(*it))
    {}

    FoamGridEntityPointer (const FoamGridEntityImp<dim-codim,dimworld>* target)
        : virtualEntity_(target)
    {}
        
        //! equality
        bool equals(const FoamGridEntityPointer<codim,GridImp>& i) const {
            return virtualEntity_.getTarget() == i.virtualEntity_.getTarget();
        }
    
        
        //! dereferencing
        Entity& dereference() const {
            return virtualEntity_;
        }

        //! ask for level of entity
        int level () const {
            return virtualEntity_.level();
        }
    
    /** \brief Throw away all temporary memory 
     * \deprecated Removed in dune-grid 2.2
     */
    void compactify() const
    {}
        
    protected:
    
        //! virtual entity
        mutable FoamGridMakeableEntity<codim,dim,GridImp> virtualEntity_;
    
        
};


} // end namespace Dune

#endif
