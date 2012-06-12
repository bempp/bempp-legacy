// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "../common/common.hpp"

#include "concrete_grid_view.hpp" // to make IDEs happy

namespace Bempp
{

template <typename DuneGridView>
void ConcreteGridView<DuneGridView>::getRawElementDataDoubleImpl(
        arma::Mat<double>& vertices,
        arma::Mat<int>& elementCorners,
        arma::Mat<char>& auxData) const
{
    getRawElementDataImpl(vertices, elementCorners, auxData);
}

template <typename DuneGridView>
void ConcreteGridView<DuneGridView>::getRawElementDataFloatImpl(
        arma::Mat<float>& vertices,
        arma::Mat<int>& elementCorners,
        arma::Mat<char>& auxData) const
{
    getRawElementDataImpl(vertices, elementCorners, auxData);
}

template <typename DuneGridView>
template <typename CoordinateType>
void ConcreteGridView<DuneGridView>::getRawElementDataImpl(
        arma::Mat<CoordinateType>& vertices,
        arma::Mat<int>& elementCorners,
        arma::Mat<char>& auxData) const
{
    typedef typename DuneGridView::Grid DuneGrid;
    typedef typename DuneGridView::IndexSet DuneIndexSet;
    const size_t dimGrid = DuneGrid::dimension;
    const size_t dimWorld = DuneGrid::dimensionworld;
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
    typedef typename DuneGrid::ctype ctype;

    const DuneIndexSet& indexSet = m_dune_gv.indexSet();

    vertices.set_size(dimWorld, indexSet.size(codimVertex));
    for (DuneVertexIterator it = m_dune_gv.template begin<codimVertex>();
         it != m_dune_gv.template end<codimVertex>(); ++it)
    {
        size_t index = indexSet.index(*it);
        const DuneVertexGeometry& geom = it->geometry();
        Dune::FieldVector<ctype, dimWorld> vertex = geom.corner(0);
        for (size_t i = 0; i < dimWorld; ++i)
            vertices(i, index) = vertex[i];
    }

    const size_t MAX_CORNER_COUNT = dimWorld == 2 ? 2 : 4;
    DuneElementMapper elementMapper(m_dune_gv.grid());
    elementCorners.set_size(MAX_CORNER_COUNT, elementMapper.size());
    for (DuneElementIterator it = m_dune_gv.template begin<codimElement>();
         it != m_dune_gv.template end<codimElement>(); ++it)
    {
        size_t index = elementMapper.map(*it);
        const Dune::GenericReferenceElement<ctype, dimGrid>& refElement =
                Dune::GenericReferenceElements<ctype, dimGrid>::general(it->type());
        const size_t cornerCount = refElement.size(codimVertex);
        assert(cornerCount <= MAX_CORNER_COUNT);
        for (size_t i = 0; i < cornerCount; ++i)
            elementCorners(i, index) = indexSet.subIndex(*it, i, codimVertex);
        for (size_t i = cornerCount; i < MAX_CORNER_COUNT; ++i)
            elementCorners(i, index) = -1;
    }

    auxData.set_size(0, elementCorners.n_cols);
}

} // namespace Bempp
