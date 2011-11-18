#ifndef DUNE_FOAMGRID_INTERSECTIONITERATORS_HH
#define DUNE_FOAMGRID_INTERSECTIONITERATORS_HH

#include <dune/foamgrid/foamgrid/foamgridintersections.hh>

/** \file
* \brief The FoamGridLeafIntersectionIterator and FoamGridLevelIntersectionIterator classes
*/

namespace Dune {

/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of an element!
*/
template<class GridImp>
class FoamGridLeafIntersectionIterator
/** \todo Inherit from the level iterator because I am too lazy to implement
    the actual leaf iterator and I don't need it yet. */
    : public FoamGridLevelIntersectionIterator<GridImp>
{
#if 0
    enum {dim=GridImp::dimension};
#endif
    
    enum {dimworld=GridImp::dimensionworld};
    
    #if 0
    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    
#endif
public:
    
    //! Constructor for a given grid entity and a given neighbor
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<2,dimworld>* center, int nb) 
        : FoamGridLevelIntersectionIterator<GridImp>(center,nb)
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLeafIntersectionIterator(const FoamGridEntityImp<2,dimworld>* center) 
        : FoamGridLevelIntersectionIterator<GridImp>(center,center->corners())
    {}

    typedef Dune::Intersection<const GridImp, Dune::FoamGridLeafIntersection> Intersection;
#if 0
    FoamGridLeafIntersectionIterator(const GridImp* identityGrid,
                                         const HostLeafIntersectionIterator& hostIterator)
        : selfLocal_(NULL), neighborLocal_(NULL), intersectionGlobal_(NULL),
          identityGrid_(identityGrid), 
          hostIterator_(hostIterator)
    {}
        
    //! equality
    bool equals(const FoamGridLeafIntersectionIterator<GridImp>& other) const {
        return hostIterator_ == other.hostIterator_;
    }

    
    //! prefix increment
    void increment() {
        ++hostIterator_;

        // Delete intersection geometry objects, if present
        if (intersectionGlobal_ != NULL) {
            delete intersectionGlobal_;
            intersectionGlobal_ = NULL;
        }
        
        if (selfLocal_ != NULL) {
            delete selfLocal_;
            selfLocal_ = NULL;
        }
        
        if (neighborLocal_ != NULL) {
            delete neighborLocal_;
            neighborLocal_ = NULL;
        }
    }
#endif
    //! \brief dereferencing
    const Intersection & dereference() const {
        return reinterpret_cast<const Intersection&>(*this);
    }
    
};




//! \todo Please doc me !
template<class GridImp>
class FoamGridLevelIntersectionIterator
{
    
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

    // Only the codim-0 entity is allowed to call the constructors
    friend class FoamGridEntity<0,dim,GridImp>;

    /** \todo Make this private once FoamGridLeafIntersectionIterator doesn't derive from this class anymore */
protected:
    //! Constructor for a given grid entity and a given neighbor
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<2,dimworld>* center, int nb) 
        : intersection_(FoamGridLevelIntersection<GridImp>(center,nb))
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    FoamGridLevelIntersectionIterator(const FoamGridEntityImp<2,dimworld>* center) 
        : intersection_(FoamGridLevelIntersection<GridImp>(center,center->corners()))
    {}

public:

    typedef Dune::Intersection<const GridImp, Dune::FoamGridLevelIntersection> Intersection;

  //! equality
  bool equals(const FoamGridLevelIntersectionIterator<GridImp>& other) const {
      return (GridImp::getRealImplementation(intersection_).center_   == GridImp::getRealImplementation(other.intersection_).center_) 
          && (GridImp::getRealImplementation(intersection_).neighbor_ == GridImp::getRealImplementation(other.intersection_).neighbor_);
  }

    //! prefix increment
    void increment() {
        GridImp::getRealImplementation(intersection_).neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
        return reinterpret_cast<const Intersection&>(intersection_);
    }

#if 0
    FoamGridElement* target() const {
        const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

        if (!isValid)
            return center_;
        else if (neighbor_==0) 
            return center_->pred_;
        else 
            return center_->succ_;

    }
#endif

private:
  //**********************************************************
  //  private data
  //**********************************************************

    /** \brief The actual intersection
    */
    mutable MakeableInterfaceObject<Intersection> intersection_;

};


}  // namespace Dune

#endif
