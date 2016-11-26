#ifndef p1_entity_pointer_hpp
#define p1_entity_pointer_hpp

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>
#include <memory>


namespace BemppGrid {

    class P1Grid;

    template<int, int, class> class P1EntityImp;
    template<int, Dune::PartitionIteratorType, class> class P1LevelIteratorImp;

    template <int codim, class>
    class P1EntityPointerImp {

        friend class P1LevelIteratorImp<codim, Dune::All_Partition, P1Grid >;

        public:

            enum {codimension = codim};
            typedef Dune::Entity<codim, 2, P1Grid, P1EntityImp> Entity;

            P1EntityPointerImp(const P1EntityImp<codim, 2, P1Grid>& entity) :
                m_entity(entity), m_duneEntity(new Entity(m_entity)) {}

            P1EntityPointerImp(const P1EntityPointerImp<codim, P1Grid>& other) :
                m_entity(other.m_entity), m_duneEntity(new Entity(m_entity)) {}

            virtual ~P1EntityPointerImp() {}

            Entity& dereference() const {

                return *m_duneEntity;

            }

            int level() const {

                return m_entity.level();

            }

            bool equals(const P1EntityPointerImp<codim, P1Grid>& other) const {

                return m_entity.equals(other.m_entity);

            }



        private:

            void setEntity(const P1EntityImp<codim, 2, P1Grid>& entity) const {

                    m_entity = entity;

                    // Need pointer as Entity has no assignment operator
                    m_duneEntity.reset(new Entity(m_entity));

            } 

            mutable P1EntityImp<codim, 2, P1Grid> m_entity;
            mutable std::unique_ptr<Entity> m_duneEntity;


        


    };



}




#endif
