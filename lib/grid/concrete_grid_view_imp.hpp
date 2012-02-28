#include "concrete_grid_view.hpp" // to make IDEs happy

namespace Bempp
{

template <typename DuneGridView>
void ConcreteGridView<DuneGridView>::getRawElementData(
        arma::Mat<ctype>& vertices,
        arma::Mat<int>& elementCorners,
        arma::Mat<char>& auxData) const
{
    typedef typename DuneGridView::Grid DuneGrid;
    typedef typename DuneGridView::IndexSet DuneIndexSet;
    const int dimGrid = DuneGrid::dimension;
    const int dimWorld = DuneGrid::dimensionworld;
    const int codimVertex = dimGrid;
    const int codimElement = 0;
    typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<DuneGrid,
            Dune::MCMGElementLayout> DuneElementMapper;
    typedef typename DuneGridView::template Codim<codimVertex>::Iterator
            DuneVertexIterator;
    typedef typename DuneGridView::template Codim<codimElement>::Iterator
            DuneElementIterator;
    typedef typename DuneGridView::template Codim<codimVertex>::Geometry
            DuneVertexGeometry;
    typedef typename DuneGridView::template Codim<codimElement>::Geometry
            DuneElementGeometry;

    const DuneIndexSet& indexSet = m_dune_gv.indexSet();

    vertices.set_size(dimWorld, indexSet.size(codimVertex));
    for (DuneVertexIterator it = m_dune_gv.template begin<codimVertex>();
         it != m_dune_gv.template end<codimVertex>(); ++it)
    {
        int index = indexSet.index(*it);
        const DuneVertexGeometry& geom = it->geometry();
        Dune::FieldVector<ctype, dimWorld> vertex = geom.corner(0);
        for (int i = 0; i < dimWorld; ++i)
            vertices(i, index) = vertex[i];
    }

    const int MAX_CORNER_COUNT = dimWorld == 2 ? 2 : 4;
    DuneElementMapper elementMapper(m_dune_gv.grid());
    elementCorners.set_size(MAX_CORNER_COUNT, elementMapper.size());
    for (DuneElementIterator it = m_dune_gv.template begin<codimElement>();
         it != m_dune_gv.template end<codimElement>(); ++it)
    {
        int index = elementMapper.map(*it);
        const Dune::GenericReferenceElement<ctype, dimGrid>& refElement =
                Dune::GenericReferenceElements<ctype, dimGrid>::general(it->type());
        const int cornerCount = refElement.size(codimVertex);
        assert(cornerCount <= MAX_CORNER_COUNT);
        for (int i = 0; i < cornerCount; ++i)
            elementCorners(i, index) = indexSet.subIndex(*it, i, codimVertex);
        for (int i = cornerCount; i < MAX_CORNER_COUNT; ++i)
            elementCorners(i, index) = -1;
    }

    auxData.set_size(0, elementCorners.n_cols);
}

} // namespace Bempp
