#ifndef bempp_p1_entity_seed_hpp
#define bempp_p1_entity_seed_hpp

namespace BemppGrid {

    class P1Grid;

    template<int, int, class> class EntityImp;

    template <int codim>
    class P1EntitySeedImp {

        public:

            enum {codimension = codim};

            P1EntitySeedImp(EntityImp<codim, 2, P1Grid> const * target ) : 
                m_target(target) {}

            bool isValid() const {

                return m_target != nullptr;

            }

        private:

            EntityImp<codim, 2, P1Grid> const * m_target;

    };

}

#endif 
