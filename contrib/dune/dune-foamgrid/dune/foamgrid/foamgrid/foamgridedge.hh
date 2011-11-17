#ifndef DUNE_FOAMGRID_EDGE_HH
#define DUNE_FOAMGRID_EDGE_HH

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>

namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<1,dimworld>
        : public FoamGridEntityBase
    {
    public:

        FoamGridEntityImp(const FoamGridEntityImp<0,dimworld>* v0, 
                          const FoamGridEntityImp<0,dimworld>* v1, 
                          int level, unsigned int id) 
            : FoamGridEntityBase(level,id)
        {
            vertex_[0] = v0;
            vertex_[1] = v1;
        }

        /** \todo Implement this method! */
        bool isLeaf() const {
            return true;
        }

        unsigned int boundarySegmentIndex() const {
            return boundaryId_;
        }

        GeometryType type() const {
            return GeometryType(1);
        }

        /** \brief Number of corners (==2) */
        unsigned int corners() const {
            return 2;
        }

        FieldVector<double, dimworld> corner(int i) const {
            return vertex_[i]->pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }
        
        // this is a helper method which only makes sense as long as
        // edges border at most two elements.
        const FoamGridEntityImp<2,dimworld>* otherElement(const FoamGridEntityImp<2,dimworld>* element) {
            assert(elements_.size()==2);
            // Return the 'other' element on the current edge
            return (elements_[0]==element) ? elements_[1] : elements_[0];
        }


        std::vector<const FoamGridEntityImp<2,dimworld>*> elements_;

        const FoamGridEntityImp<0,dimworld>* vertex_[2];

        /** \brief The boundary id.  Only used if this edge is a boundary edge */
        unsigned int boundaryId_;
        
    };

}

#endif
