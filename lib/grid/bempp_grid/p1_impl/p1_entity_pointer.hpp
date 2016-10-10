#ifndef p1_entity_pointer_hpp
#define p1_entity_pointer_hpp

#include <dune/grid/common/entity.hh>

namespace BemppGrid {

    class P1Grid;

    template<int, int, class> class P1EntityImp;

    template <int codim>
    class P1EntityPointerImp {

        public:

            enum {codimension = codim};
            typedef Dune::Entity<codim, 2, P1Grid, P1EntityImp> Entity;

            P1EntityPointerImp(const P1EntityImp<codim, 2, P1Grid>& entity) :
                m_entity(entity), m_duneEntity(entity) {}

            P1EntityPointerImp(const P1EntityPointerImp<codim>& other) :
                m_entity(other.m_entity), m_duneEntity(other.m_entity) {}

            Entity& dereference() const {

                return m_duneEntity;

            }

            int level() const {

                return m_entity.level();

            }

            bool equals(const P1EntityPointerImp<codim>& other) {

                return m_entity.equals(other.m_entity);

            }


        private:

            P1EntityImp<codim, 2, P1Grid> m_entity;
            mutable Entity m_duneEntity;


        


    };



}




#endif
