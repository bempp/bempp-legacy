#ifndef bempp_p1_entity_seed_hpp
#define bempp_p1_entity_seed_hpp

namespace BemppGrid {

    class P1Grid;

    template<int, int, class> class P1EntityImp;

    template <int codim, class>
    class P1EntitySeedImp {

        public:

            enum {codimension = codim};

            P1EntitySeedImp(int level, std::size_t index) : 
                m_level(level), m_index(index) {}

            bool isValid() const {

                return true;

            }

        private:

            int m_level;
            std::size_t m_index;


    };

}

#endif 
