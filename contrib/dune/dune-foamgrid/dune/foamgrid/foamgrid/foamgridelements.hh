// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=8 sw=4 et sts=4:
#ifndef DUNE_FOAMGRID_ELEMENTS_HH
#define DUNE_FOAMGRID_ELEMENTS_HH

#include <dune/foamgrid/foamgrid/foamgridvertex.hh>
#include <dune/foamgrid/foamgrid/foamgridedge.hh>
#include <dune/common/nullptr.hh>
namespace Dune {

    template <int dimworld>
    class FoamGridEntityImp<2,dimworld>
        : public FoamGridEntityBase
    {
    public:

        /** \brief The different ways to mark an element for grid changes */
        enum MarkState { DO_NOTHING , COARSEN , REFINE, IS_COARSENED };

        FoamGridEntityImp(int level, unsigned int id) 
            : FoamGridEntityBase(level,id),
              nSons_(0), refinementIndex_(-1),
              markState_(DO_NOTHING), isNew_(false)
        {
          sons_[0]= sons_[1] = sons_[2] = sons_[3] = nullptr;
          father_ = nullptr;
        }

        bool hasFather() const
        {
            return father_!=nullptr;
        }
        
        bool mightVanish() const
        {
            return markState_==COARSEN;
        }
        
        bool isLeaf() const {
            // The sons are either all nullptr or all != nullptr
            return sons_[0] == nullptr;
        }

        bool isNew() const 
        {
            return isNew_;
        }
        
        unsigned int corners() const {
            return 3;
        }

        GeometryType type() const {
            return GeometryType(GeometryType::simplex, 2);
        }

        /** \todo Implement me! */
        unsigned int nSons() const {
            return nSons_;
        }
      
        unsigned int nSons_;

        /** \brief Compute local cordinates from global ones.
         * \param coord The global coordinates.
         * \return The corresponding local coordinates within the element.
         */
        FieldVector<double,2> globalToLocal(const FieldVector<double,dimworld>& coord) const
        {
            // If we set up the overdetermined system matrix we have
            // A[i][0]=vertex_[1].pos_[i]-vertex_[0].pos_[i];
            // A[i][1]=vertex_[2].pos_[i]-vertex_[0].pos_[i];
            // t[i]=coord[i]-vertex_[0].pos_[i];
            //
            // to determine the local coordinates we solve
            // A'A x= A' t
            //
            
            FieldMatrix<double,2,2> mat; // A'A
            // mat_{ij}=\sum_k A_{ki}A_{kj}
            mat=0;
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[0][0]+=(vertex_[1]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[1]->pos_[i]-vertex_[0]->pos_[i]);
            }
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[1][0]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[1]->pos_[i]-vertex_[0]->pos_[i]);
            }
            mat[0][1]=mat[1][0];
            for(std::size_t i=0; i <dimworld; ++i)
            {
                mat[1][1]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(vertex_[2]->pos_[i]-vertex_[0]->pos_[i]);
            }
            
            FieldVector<double, 2> b, x;
            b=0;
            for(std::size_t i=0; i <dimworld; ++i)
            {
                b[0]+=(vertex_[1]->pos_[i]-vertex_[0]->pos_[i])*(coord[i]-vertex_[0]->pos_[i]);
                b[1]+=(vertex_[2]->pos_[i]-vertex_[0]->pos_[i])*(coord[i]-vertex_[0]->pos_[i]);
            }
            mat.solve(x, b);
#ifndef NDEBUG
            FieldVector<double,dimworld> test(vertex_[0]->pos_);
            test.axpy(x[0], vertex_[1]->pos_);
            test.axpy(-x[0], vertex_[0]->pos_);
            test.axpy(x[1], vertex_[2]->pos_);
            test.axpy(-x[1], vertex_[0]->pos_);
            assert((test-coord).two_norm()< std::numeric_limits<double>::epsilon()*8);
#endif
            return x;
        }
        
        /**
         * \brief index of the refined element in the father
         *
         * For red refinement this is either the index of corner,
         * that is also a corner in the father element, within the father
         * or 3 if no corner is also a corner in the father.
         */
        int refinementIndex_;

        array<FoamGridEntityImp<2,dimworld>*, 4> sons_;

        FoamGridEntityImp<2,dimworld>* father_;

        array<FoamGridEntityImp<1,dimworld>*, 3> edges_;

        FoamGridEntityImp<0,dimworld>* vertex_[3];
        
        /** \brief Stores requests for refinement and coarsening */
        MarkState markState_;
        
        /** \brief This flag is set by adapt() if this element has been newly created. */
        bool isNew_;
        
    };

}

#endif
