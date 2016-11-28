#ifndef bempp_grid_triangle_imp_entity_seed_hpp
#define bempp_grid_triangle_imp_entity_seed_hpp

namespace BemppGrid {

    class TriangleGrid;

    template<int, int, class> class EntityImp;

    template <int codim, class>
    class EntitySeedImp {

        public:

            enum {codimension = codim};

            EntitySeedImp(int level, std::size_t index) : 
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
