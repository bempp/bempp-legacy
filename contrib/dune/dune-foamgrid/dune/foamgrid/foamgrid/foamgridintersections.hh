// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_INTERSECTIONS_HH
#define DUNE_FOAMGRID_INTERSECTIONS_HH

/** \file
* \brief The FoamGridLeafIntersection and FoamGridLevelIntersection classes
*/

#include <dune/grid/common/intersection.hh>

#include <dune/foamgrid/foamgrid/foamgridintersectioniterators.hh>
#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridgeometry.hh>
#include <dune/foamgrid/foamgrid/foamgridentitypointer.hh>

namespace Dune {

template <class GridImp>
class FoamGridLevelIntersectionIterator;


//! \brief Base class of all intersections within FoamGrid
//!
//! encapsulates common functionality of level and leaf intersections.
template<class GridImp>
class FoamGridIntersection
{
        
        // The type used to store coordinates
        typedef typename GridImp::ctype ctype;

    friend class FoamGridLevelIntersectionIterator<GridImp>;
    friend class FoamGridLeafIntersectionIterator<GridImp>;
    
    public:

    
        enum {dim=GridImp::dimension};
    
        enum {dimworld=GridImp::dimensionworld};

        typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
        typedef typename GridImp::template Codim<1>::Geometry Geometry;
        typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
        typedef typename GridImp::template Codim<0>::Entity Entity;
        typedef Dune::Intersection<const GridImp, Dune::FoamGridLevelIntersectionIterator> Intersection;

    /**
     * \brief Initalizes an intersection.
     *
     * After initialization this object always represents the first intersection 
     * related to an edge.
     * \param edge The index of the edge this intersection lives on.
     */
    FoamGridIntersection(const FoamGridEntityImp<2,dimworld>* center,
                         int edge)
        : center_(center), edgeIndex_(edge)
    {}
    /*
    void increment(){
        ++neighbor_;
        while(neighbor_<(*edges_)[edgeIndex_]->elements_.size() &&
              (center_==center_->edges_[edgeIndex_]->elements_[neighbor_]
               ||center_->level()!=center_->edges_[edgeIndex_]->elements_[neighbor_]->level()))
            ++neighbor_;
        
        if(neighbor_<center_->edges_[edgeIndex_]->elements_.size())
            return;
        neighbor_=0;
        for(++edgeIndex_;edgeIndex_!=center_->corners();++edgeIndex_){ // Not an end iterator
            while(center_->edges_[edgeIndex_]->elements_.size()>1 // not on boundary
                  && (center_==center_->edges_[edgeIndex_]->elements_[neighbor_]
                      ||center_->level()!=center_->edges_[edgeIndex_]->elements_[neighbor_]->level()))
            {
                // Move index to point to the first real neighbor 
                ++neighbor_;
            }
            if(neighbor_<center_->edges_[edgeIndex_]->elements_.size())
            {
                // Found the first intersection
                break;
            }
            neighbor_=0;
        }
    }
    */
    
        //! return EntityPointer to the Entity on the inside of this intersection
        //! (that is the Entity where we started this Iterator)
        EntityPointer inside() const {
            return FoamGridEntityPointer<0,GridImp> (center_);
        }

        
        //! return EntityPointer to the Entity on the outside of this intersection
        //! (that is the neighboring Entity)
        EntityPointer outside() const {
            // Return the 'other' element on the current edge
            return FoamGridEntityPointer<0,GridImp> ((*neighbor_));
        }
        
        
    /** \brief return true if intersection is with boundary.
    */
    bool boundary () const {
        return center_->edges_[edgeIndex_]->elements_.size()==1;
    }
        
        
    //! return information about the Boundary
    int boundaryId () const DUNE_DEPRECATED {
        return center_->edges_[edgeIndex_]->boundarySegmentIndex();
    }
        
    //! return information about the Boundary
    int boundarySegmentIndex () const {
        return center_->edges_[edgeIndex_]->boundarySegmentIndex();
    }
        
    //! Return true if this is a conforming intersection
    bool conforming () const {
        // FoamGrid level intersections are always conforming
        return true;
    }
        
    //! Geometry type of an intersection
    GeometryType type () const {
        return GeometryType(GeometryType::simplex, dim-1);
    }


        //! intersection of codimension 1 of this neighbor with element where
        //! iteration started.
        //! Here returned element is in LOCAL coordinates of the element
        //! where iteration started.
        LocalGeometry geometryInInside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            // Get two vertices of the intersection
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            coordinates[0] = refElement.position(refElement.subEntity(edgeIndex_, 1, 0, dim),dim);
            coordinates[1] = refElement.position(refElement.subEntity(edgeIndex_, 1, 1, dim),dim);

            return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(type(), coordinates));
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in LOCAL coordinates of neighbor
        LocalGeometry geometryInOutside () const {

            std::vector<FieldVector<double, dim> > coordinates(2);

            // Get two vertices of the intersection
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general((*neighbor_)->type());

            int idxInOutside = indexInOutside();
            
            coordinates[0] = refElement.position(refElement.subEntity(idxInOutside, 1, 0, dim),dim);
            coordinates[1] = refElement.position(refElement.subEntity(idxInOutside, 1, 1, dim),dim);

            return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(type(), coordinates));
        }
        
        //! intersection of codimension 1 of this neighbor with element where iteration started.
        //! Here returned element is in GLOBAL coordinates of the element where iteration started.
        Geometry geometry () const {

            std::vector<FieldVector<double, dimworld> > coordinates(2);

            // Get two vertices of the intersection
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            coordinates[0] = center_->vertex_[refElement.subEntity(edgeIndex_, 1, 0, dim)]->pos_;
            coordinates[1] = center_->vertex_[refElement.subEntity(edgeIndex_, 1, 1, dim)]->pos_;

            return Geometry(FoamGridGeometry<dim-1, dimworld, GridImp>(type(), coordinates));
        }
        
        
        //! local number of codim 1 entity in self where intersection is contained in
        int indexInInside () const {
            return edgeIndex_;
        }
        
        virtual int indexInOutside() const=0;
    
        //! return outer normal
        FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
            // The intersection normal is a vector that is orthogonal to the element normal
            // and to the intersection itself.
            
            // only works for triangles
            assert(center_->type().isTriangle());

            // Compute vertices
            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            // edge vertices, oriented
            int v0 = std::min(refElement.subEntity(edgeIndex_, 1, 0, dim), refElement.subEntity(edgeIndex_, 1, 1, dim));
            int v1 = std::max(refElement.subEntity(edgeIndex_, 1, 0, dim), refElement.subEntity(edgeIndex_, 1, 1, dim));

            // opposite vertex
            int v2 = (v1+1)%3;

            // Compute oriented edge
            FieldVector<ctype, dimworld> edge = center_->vertex_[v1]->pos_ - center_->vertex_[v0]->pos_;

            // compute triangle edge normal
            FieldVector<ctype, dimworld> scaledEdge = edge;
            edge *= edge*(center_->vertex_[v2]->pos_ - center_->vertex_[v0]->pos_);
            FieldVector<ctype, dimworld> normal = center_->vertex_[v2]->pos_ - center_->vertex_[v0]->pos_;
            normal -= scaledEdge;
            normal *= -1;
            return normal;
        }

        //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {

            const Dune::GenericReferenceElement<double,dim>& refElement
                = Dune::GenericReferenceElements<double, dim>::general(center_->type());

            ctype edgeLength = (center_->vertex_[refElement.subEntity(edgeIndex_, 1, 0, dim)]->pos_
                                - center_->vertex_[refElement.subEntity(edgeIndex_, 1, 1, dim)]->pos_).two_norm();

            FieldVector<ctype, dimworld> normal = unitOuterNormal(local);
            normal *= edgeLength;
            return normal;
        }

        //! return unit outer normal
        FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
            FieldVector<ctype, dimworld> outerNormal = this->outerNormal(local);
            outerNormal /= outerNormal.two_norm();
            return outerNormal;
        }

        //! return unit outer normal at the intersection center
        FieldVector<ctype, dimworld> centerUnitOuterNormal () const {
            FieldVector<ctype, dimworld> outerNormal = this->outerNormal(FieldVector<ctype,1>(0.5));
            outerNormal /= outerNormal.two_norm();
            return outerNormal;
        }
    private:
 
    //! vector storing the outer normal 
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    protected:

    const FoamGridEntityImp<2,dimworld>* center_;

    /** \brief Count on which edge we are lookin' at.  */
    int edgeIndex_;
    
    /** \brief Iterator to the other neighbor of the intersection. */
    typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator neighbor_;
    
};

template<class GridImp>
class FoamGridLevelIntersection
    : public FoamGridIntersection<GridImp>
{
    friend class FoamGridLevelIntersectionIterator<GridImp>;

public:
    enum{ dimworld = GridImp::dimensionworld };
    
    FoamGridLevelIntersection(const FoamGridEntityImp<2,dimworld>* center, int edge)
                              : FoamGridIntersection<GridImp>(center, edge)
    {}

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
        assert(this->center_->edges_[this->edgeIndex_]->elements_.size()==2);
        assert(this->neighborIndex_!=this->center_->edges_[this->edgeIndex_]->elements_.size());
        
        return std::find((*this->neighbor_)->edges_.begin(), (*this->neighbor_)->edges_.end(), 
                         this->center_->edges_[this->edgeIndex_]) 
            - (*this->neighbor_)->edges_.begin();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return this->neighborIndex_!=this->center_->edges_[this->edgeIndex_]->elements_.size();
    }
    private:
    int neighborIndex_;
    
};
    
        


/** \brief Iterator over all element neighbors
* \ingroup FoamGrid
* Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
* a neighbor is an entity of codimension 0 which has a common entity of codimension 1
* These neighbors are accessed via a IntersectionIterator. This allows the implement
* non-matching meshes. The number of neighbors may be different from the number
* of an element!
*/
template<class GridImp>
class FoamGridLeafIntersection
    : public FoamGridIntersection<GridImp>
{

    friend class FoamGridLeafIntersectionIterator<GridImp>;
public:
    
    enum {dimworld=GridImp::dimensionworld};
    
    enum {dim=GridImp::dimension};

    typedef Dune::Intersection<const GridImp, Dune::FoamGridLeafIntersection> Intersection;

    typedef typename GridImp::ctype ctype;
    
    typedef typename FoamGridIntersection<GridImp>::LocalGeometry LocalGeometry;

    FoamGridLeafIntersection(const FoamGridEntityImp<2,FoamGridIntersection<GridImp>::dimworld>* center,
                             int edge)
        : FoamGridIntersection<GridImp>(center, edge)
    {}    
    
    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
        assert(this->neighbor_!=neighborEnd_);
        // Move to the father of the edge until its level is the same as
        // the level of the neighbor
        FoamGridEntityImp<1,dimworld>* edge=(this->center_->edges_[this->edgeIndex_]);
        
        while(edge->level()<(*this->neighbor_)->level())
        {
            assert(edge->father_!=nullptr);
            edge=edge->father_;
        }
        assert(edge->level()==(*this->neighbor_)->level());
        assert(edge->elements_.size()==2);
        return std::find((*this->neighbor_)->edges_.begin(), (*this->neighbor_)->edges_.end(), edge) 
            - (*this->neighbor_)->edges_.begin();
        
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const {
        
        if(edgePointer_->level()==this->center_->level())
            return FoamGridIntersection<GridImp>::geometryInInside();
        
        std::vector<FieldVector<double, dim> > coordinates(2);
        
        coordinates[0]=this->center_->globalToLocal(edgePointer_->vertex_[0]->pos_);
        coordinates[0]=this->center_->globalToLocal(edgePointer_->vertex_[1]->pos_);
        
        return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInOutside () const {
        /*
        if(edgePointer_.level()==(*this->neighbor_)->level())
            return FoamGridIntersection<GridImp>::geometryInOutside();
        */
        std::vector<FieldVector<double, dim> > coordinates(2);
        
        coordinates[0]=(*this->neighbor_)->globalToLocal(edgePointer_->vertex_[0]->pos_);
        coordinates[0]=(*this->neighbor_)->globalToLocal(edgePointer_->vertex_[1]->pos_);
        
        return LocalGeometry(FoamGridGeometry<dim-1, dim, GridImp>(this->type(), coordinates));
    }
  
     //! return outer normal multiplied by the integration element
        FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
            ctype edgeLength = (edgePointer_->vertex_[0]->pos_ - edgePointer_->vertex_[1]->pos_).two_norm();
            FieldVector<ctype, dimworld> normal = this->unitOuterNormal(local);
            normal *= edgeLength;
            return normal;
        }
    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
        return this->neighbor_!=neighborEnd_;
    }
private:
    FoamGridEntityImp<1,dimworld>* edgePointer_;
    /** \brief Iterator to the other neighbor of the intersection. */
    typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator neighborEnd_;
};


}  // namespace Dune

#endif
