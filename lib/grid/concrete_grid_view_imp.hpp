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
#ifndef bempp_concrete_grid_view_imp_hpp
#define bempp_concrete_grid_view_imp_hpp

#include "../common/acc.hpp"
#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "concrete_grid_view.hpp" // to make IDEs happy

namespace Bempp {

template <typename DuneGridView>
void ConcreteGridView<DuneGridView>::getRawElementDataDoubleImpl(
    Matrix<double> &vertices, Matrix<int> &elementCorners,
    Matrix<char> &auxData, std::vector<int> *domainIndices) const {
  getRawElementDataImpl(vertices, elementCorners, auxData, domainIndices);
}

template <typename DuneGridView>
void ConcreteGridView<DuneGridView>::getRawElementDataFloatImpl(
    Matrix<float> &vertices, Matrix<int> &elementCorners, Matrix<char> &auxData,
    std::vector<int> *domainIndices) const {
  getRawElementDataImpl(vertices, elementCorners, auxData, domainIndices);
}

template <typename DuneGridView>
template <typename CoordinateType>
void ConcreteGridView<DuneGridView>::getRawElementDataImpl(
    Matrix<CoordinateType> &vertices, Matrix<int> &elementCorners,
    Matrix<char> &auxData, std::vector<int> *domainIndices) const {
  typedef typename DuneGridView::Grid DuneGrid;
  typedef typename DuneGridView::IndexSet DuneIndexSet;
  const int dimGrid = DuneGrid::dimension;
  const int dimWorld = DuneGrid::dimensionworld;
  const int codimVertex = dimGrid;
  const int codimElement = 0;
  typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<DuneGrid,
                                                        Dune::MCMGElementLayout>
      DuneElementMapper;
  typedef typename DuneGridView::template Codim<codimVertex>::Iterator
      DuneVertexIterator;
  typedef typename DuneGridView::template Codim<codimElement>::Iterator
      DuneElementIterator;
  typedef typename DuneGridView::template Codim<codimVertex>::Geometry
      DuneVertexGeometry;
  typedef typename DuneGridView::template Codim<codimElement>::Geometry
      DuneElementGeometry;
  typedef typename DuneGrid::ctype ctype;

  const DuneIndexSet &indexSet = m_dune_gv.indexSet();

  vertices.resize(dimWorld, indexSet.size(codimVertex));
  for (DuneVertexIterator it = m_dune_gv.template begin<codimVertex>();
       it != m_dune_gv.template end<codimVertex>(); ++it) {
    size_t index = indexSet.index(*it);
    const DuneVertexGeometry &geom = it->geometry();
    Dune::FieldVector<ctype, dimWorld> vertex = geom.corner(0);
    for (int i = 0; i < dimWorld; ++i)
      vertices(i, index) = vertex[i];
  }

  const int MAX_CORNER_COUNT = dimWorld == 2 ? 2 : 4;
  DuneElementMapper elementMapper(m_dune_gv.grid());
  const int elementCount = elementMapper.size();
  elementCorners.resize(MAX_CORNER_COUNT, elementCount);
  for (DuneElementIterator it = m_dune_gv.template begin<codimElement>();
       it != m_dune_gv.template end<codimElement>(); ++it) {
    size_t index = indexSet.index(*it);
    const Dune::ReferenceElement<ctype, dimGrid> &refElement =
        Dune::ReferenceElements<ctype, dimGrid>::general(it->type());
    const int cornerCount = refElement.size(codimVertex);
    assert(cornerCount <= MAX_CORNER_COUNT);
    for (int i = 0; i < cornerCount; ++i)
      elementCorners(i, index) = indexSet.subIndex(*it, i, codimVertex);
    for (int i = cornerCount; i < MAX_CORNER_COUNT; ++i)
      elementCorners(i, index) = -1;
  }

  auxData.resize(0, elementCorners.cols());

  if (domainIndices) {
    // Somewhat inelegant: we perform a second iteration over elements,
    // this time using the BEM++ interface to Dune.
    domainIndices->resize(elementCount);
    std::unique_ptr<EntityIterator<0>> it = this->entityIterator<0>();
    const IndexSet &indexSet = this->indexSet();
    while (!it->finished()) {
      const Entity<0> &entity = it->entity();
      const int index = indexSet.entityIndex(entity);
      const int domain = entity.domain();
      acc(*domainIndices, index) = domain;
      it->next();
    }
  }
}

} // namespace Bempp
#endif
