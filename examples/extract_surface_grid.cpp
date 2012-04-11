// Copyright (C) 2011 by the BEM++ Authors
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

#include <armadillo>
#include <deque>
#include <iostream>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "grid/concrete_grid.hpp"
#include "grid/dune.hpp"
#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"

using namespace Bempp;

template <typename DuneGrid>
void doStuff()
{
    const char MESH_FNAME[] = "full-sphere.gmsh";
    const int ANY_PHYSICAL_ENTITY = -1;
    const int requestedPhysicalEntity = 1;

    // const int worldDim = DuneGrid::worlddimension;
    const int gridDim = DuneGrid::dimension;
    const int vertexCodim = gridDim;
    const int faceCodim = 1;

    GridParameters params;
    params.topology = GridParameters::TETRAHEDRAL;

    // probably this should be renamed to boundarySegmentIndex2...
    std::vector<int> boundaryId2PhysicalEntity, elementIndex2PhysicalEntity;

    std::auto_ptr<Grid> grid(
                GridFactory::importGmshGrid(
                    params, std::string(MESH_FNAME),
                    boundaryId2PhysicalEntity, elementIndex2PhysicalEntity,
                    true, // verbose
                    false)); // insertBoundarySegments

    const DuneGrid& duneGrid =
            dynamic_cast<const ConcreteGrid<DuneGrid>&>(*grid).duneGrid();

    typedef typename DuneGrid::LeafGridView DuneGridView;
    DuneGridView view = duneGrid.leafView();

    typedef typename DuneGridView::IndexSet DuneIndexSet;
    typedef typename DuneIndexSet::IndexType IndexType;
    const DuneIndexSet& indexSet = view.indexSet();

    std::set<IndexType> surfaceVertexIndices;
    std::deque<std::pair<GeometryType, std::vector<IndexType> > > surfaceElements;

    typedef typename DuneGridView::template Codim<0>::Iterator DuneElementIterator;
    for (DuneElementIterator eit = view.template begin<0>();
         eit != view.template end<0>(); ++eit)
    {
        typedef typename DuneGridView::template Codim<0>::Entity DuneElement;
        const DuneElement& element = *eit;

        typedef typename DuneGridView::IntersectionIterator DuneIntersectionIterator;
        for (DuneIntersectionIterator iit = view.ibegin(element);
             iit != view.iend(element); ++iit)
        {
            typedef typename DuneGridView::Intersection DuneIntersection;
            const DuneIntersection& intersection = *iit;
            if (!intersection.boundary())
                continue;
            int physicalEntity =
                    boundaryId2PhysicalEntity[intersection.boundarySegmentIndex()];
            if (requestedPhysicalEntity != ANY_PHYSICAL_ENTITY &&
                    requestedPhysicalEntity != physicalEntity)
                continue;

            // Get the indices of the vertices of the face which constitutes the intersection.
            // This is a slightly complicated.

            int faceLocalIndex = intersection.indexInInside();

            typedef Dune::GenericReferenceElement<ctype, gridDim> GRE;
            typedef Dune::GenericReferenceElements<ctype, gridDim> GREs;
            const GRE& refElement = GREs::general(element.type());
            GeometryType faceType = refElement.type(faceLocalIndex, faceCodim);
            int vertexCount =
                    refElement.size(faceLocalIndex, faceCodim, vertexCodim);

            // Global indices of the vertices making up the face
            std::vector<IndexType> vertexIndices(vertexCount);
            for (int i = 0; i < vertexCount; ++i)
            {
                int vertexLocalIndex = refElement.subEntity(
                            faceLocalIndex, faceCodim, i, vertexCodim);
                IndexType vertexGlobalIndex =
                        indexSet.subIndex(element, vertexLocalIndex, vertexCodim);
                surfaceVertexIndices.insert(vertexGlobalIndex);
                vertexIndices[i] = vertexGlobalIndex;
            }
            // Note: one might try to optimise away this vector copy
            surfaceElements.push_back(std::make_pair(faceType, vertexIndices));

            // UNFINISHED -- TO BE CONTINUED
        }
    }
}

int main()
{
    doStuff<Default3dIn3dDuneGrid>();
}

