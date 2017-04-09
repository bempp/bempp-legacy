#ifndef bempp_grid_imp_triangle_factory_hpp
#define bempp_grid_imp_triangle_factory_hpp

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>
#include "bempp_grid_data_container.hpp"
#include <memory>

namespace Dune {

template <>
class GridFactory<Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>>
    : public Dune::GridFactoryInterface<
          Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily>> {

public:
  enum { dimworld = 3 };
  typedef double ctype;

  GridFactory()
      : m_nodes(new std::vector<Dune::FieldVector<double, 3>>),
        m_elements(new std::vector<std::vector<unsigned int>>) {}

  virtual void
  insertVertex(const Dune::FieldVector<ctype, dimworld> &pos) override {

    m_nodes->push_back(pos);
  }

  virtual void
  insertElement(const Dune::GeometryType &type,
                const std::vector<unsigned int> &vertices) override {

    m_elements->push_back(vertices);
  }

  virtual Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily> *
  createGrid() override {

    auto container = BemppGrid::shared_ptr<BemppGrid::DataContainer>(
        new BemppGrid::DataContainer());
    container->init(m_nodes, m_elements);
    m_grid = new BemppGrid::TriangleGrid(container);
    return m_grid;
  }

  virtual unsigned int insertionIndex(
      const typename BemppGrid::TriangleGridFamily::Traits::Codim<0>::Entity
          &entity) const override {

    assert(entity.level() == 0);
    return m_grid->levelIndexSet(0).index<0>(entity);
  }

  virtual unsigned int insertionIndex(
      const typename BemppGrid::TriangleGridFamily::Traits::Codim<2>::Entity
          &entity) const override {

    assert(entity.level() == 0);
    return m_grid->levelIndexSet(0).index<2>(entity);
  }

  virtual void
  insertBoundarySegment(const std::vector<unsigned int> &vertices) override {

    throw std::runtime_error(
        "GridFactory::insertBoundarySegment(): Not implemented.");
  }

  virtual void
  insertBoundarySegment(const std::vector<unsigned int> &vertices,
                        const std::shared_ptr<BoundarySegment<dimension, dimworld>>
                            &boundarySegment) override {

    std::runtime_error(
        "Gridfactory::insertBoundarySegment(): not implemented.");
  }

private:
  BemppGrid::shared_ptr<std::vector<Dune::FieldVector<double, 3>>> m_nodes;
  BemppGrid::shared_ptr<std::vector<std::vector<unsigned int>>> m_elements;
  Dune::Grid<2, 3, double, BemppGrid::TriangleGridFamily> *m_grid;
};
}

#endif
