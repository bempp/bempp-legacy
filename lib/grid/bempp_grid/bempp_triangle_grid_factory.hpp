#ifndef bempp_triangle_grid_factory_hpp
#define bempp_triangle_grid_factory_hpp

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>
#include "./triangle_imp/bempp_grid_data_container.hpp"

namespace Dune {

    template <>
    class GridFactory<Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>> : 
        public Dune::GridFactoryInterface<Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>> {


            public:

                enum {dimworld = 3};
                typedef double ctype;

                GridFactory() :
                    m_nodes(new std::vector<Dune::FieldVector<double, 3>>),
                    m_elements(new std::vector<std::vector<unsigned int>>) {}

                void insertVertex(const Dune::FieldVector<ctype, dimworld>& pos) {

                    m_nodes->push_back(pos);

                }

                void insertElement(const Dune::GeometryType& type, const std::vector<unsigned int>& vertices) {


                    m_elements->push_back(vertices);

                }

                Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>* createGrid() {


                    auto container = BemppGrid::shared_ptr<BemppGrid::DataContainer>(new BemppGrid::DataContainer());
                    container->init(m_nodes, m_elements);
                    m_grid = new BemppGrid::TriangleGrid(container);
                    return m_grid;
                }

                unsigned int insertionIndex(const typename BemppGrid::TriangleGridFamily::Traits::Codim<0>::Entity& entity) const {

                    assert(entity.level() == 0);
                    return m_grid->levelIndexSet(0).index<0>(entity);

                }

                unsigned int insertionIndex(const typename BemppGrid::TriangleGridFamily::Traits::Codim<2>::Entity& entity) const {

                    assert(entity.level() == 0);
                    return m_grid->levelIndexSet(0).index<2>(entity);

                }

                void insertBoundarySegment(const std::vector<unsigned int>& vertices) {

                    throw std::runtime_error("GridFactory::insertBoundarySegment(): Not implemented.");

                }

            private:

                BemppGrid::shared_ptr<std::vector<Dune::FieldVector<double, 3>>> m_nodes;
                BemppGrid::shared_ptr<std::vector<std::vector<unsigned int>>> m_elements;
                Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>* m_grid;

        };


}


#endif
