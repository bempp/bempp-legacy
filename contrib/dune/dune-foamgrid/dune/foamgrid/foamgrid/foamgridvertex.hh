// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_VERTEX_HH
#define DUNE_FOAMGRID_VERTEX_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>


namespace Dune {

    /** \brief Base class for FoamGrid entity implementation classes */
    class FoamGridEntityBase
    {
    public:
        FoamGridEntityBase(int level, unsigned int id)
            : level_(level), id_(id), willVanish_()
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
        //! \brief Whether this entity will vanish due to coarsening.
        bool willVanish_;
    };

    /**
     * \brief The actual entity implementation
     *
     * \tparam dim The dimension of this entity
     * \tparam dimworld The world diemnsion
     */
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
              pos_(pos), son_(nullptr) 
        {}
        
        //private: 
        bool isLeaf() const {
            return son_==nullptr;
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
