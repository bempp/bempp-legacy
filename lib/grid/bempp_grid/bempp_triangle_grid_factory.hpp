#ifndef bempp_triangle_grid_factory_hpp
#define bempp_triangle_grid_factory_hpp

#include <dune/grid/common/gridfactory.hh>
#include "./triangle_imp/bempp_grid_data_container.hpp"

namespace BemppGrid {

    template <>
    Dune::GridFactory<Dune::Grid<2, 3, double, TriangleGridFamily>> : 
        public Dune::GridFactoryInterface<Dune::Grid<2, 3, double, TriangleGridFamily>> {


            public:

                enum {dimworld = 3};
                typedef double ctype;

                GridFactory() :
                    m_nodes(new std::vector<Dune::FieldVector<double, 3>>),
                    m_elements(new std::vector<std::vector<unsigned int>>> m_elements) {}

                void insertVertex(const Dune::FieldVector<ctype, dimworld>& pos) {

                    m_nodes->push_back(pos);

                }

                void insertElement(const Dune::GeometryType& type, const std::vector<unsigned int>& vertices) {


                    m_elements->push_back(vertices);

                }

                Dune::Grid<2, 3, double, TriangleGridFamily>* GridFactory<Dune::Grid<2, 3, double, TriangleGridFamily>>::createGrid() {


                    auto container = shared_ptr<DataContainer>(new DataContainer());
                    container->init(m_nodes, m_elements);
                    m_grid = new TriangleGrid(container);
                    return m_grid;
                }

                unsigned int insertionIndex(const typename TriangleGridFamily::Traits::Codim<0>::Entity& entity) const {

                    assert(entity.level() == 0);
                    return m_grid->levelIndexSet(0).index<0>(entity);

                }

                unsigned int insertionIndex(const typename TriangleGridFamily::Traits::Codim<2>::Entity& entity) const {

                    assert(entity.level() == 0);
                    return m_grid->levelIndexSet(0).index<2>(entity);

                }

            private:

                shared_ptr<std::vector<Dune::FieldVector<double, 3>>> m_nodes;
                shared_ptr<std::vector<std::vector<unsigned int>>> m_elements;
                GridType* m_grid;

        };


}


#endif
