#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

namespace Dune {

    /** \brief Base class for FoamGrid entity implementation classes */
    class FoamGridEntityBase
    {
    public:
        FoamGridEntityBase(int level, unsigned int id)
            : level_(level), id_(id)
        {}

        unsigned int level() const {
            return level_;
        }

        //! level
        int level_;
        
        //! entity number 
        unsigned int levelIndex_;
        
        unsigned int leafIndex_;
        
        unsigned int id_;
        
    };

    template <int dim, int dimworld>
    class FoamGridEntityImp
    {};

    template <int dimworld>
    class FoamGridEntityImp<0,dimworld>
        : public FoamGridEntityBase
    {
    public:
        
        FoamGridEntityImp(int level, const FieldVector<double, dimworld>& pos, unsigned int id) 
            : FoamGridEntityBase(level, id),
              pos_(pos), son_(NULL) 
        {}
        
        //private: 
        bool isLeaf() const {
            return son_==NULL;
        }

        GeometryType type() const {
            return GeometryType(0);
        }

        /** \brief Number of corners (==1) */
        unsigned int corners() const {
            return 1;
        }

        FieldVector<double, dimworld> corner(int i) const {
            return pos_;
        }

        PartitionType partitionType() const {
            return InteriorEntity;
        }
        
        FieldVector<double, dimworld> pos_;
        
        //! Son vertex on the next finer grid
        FoamGridEntityImp<0,dimworld>* son_;
        
    };

}

#endif
