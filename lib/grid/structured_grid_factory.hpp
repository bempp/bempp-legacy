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

#ifndef bempp_structured_grid_factory_hpp
#define bempp_structured_grid_factory_hpp

#include "../common/common.hpp"

/** \file
    \brief A class to construct structured cube and simplex grids using the grid factory
*/

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <memory>

#include <dune/common/array.hh>
#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>

namespace Dune
{

/** \brief Construct structured cube and simplex grids in unstructured grid managers

\note Copied from Dune and slightly modified:
\li Grids embedded in higher-dimensional spaces are supported. In this case the last (dimworld - dim)
  coordinates of all vertices, where dim is the grid dimension and dimworld the world dimension,
  are set to zero.
\li Automatic pointers (std::auto_ptr) instead of shared pointers (Dune::shared_ptr) are returned.
*/
template <class GridType>
class BemppStructuredGridFactory
{
    typedef typename GridType::ctype ctype;

    static const int dim = GridType::dimension;

    static const int dimworld = GridType::dimensionworld;

    /** \brief dim-dimensional multi-index.  The range for each component can be set individually
    */
    class MultiIndex
        : public array<int,dim>
    {

        // The range of each component
        array<int,dim> limits_;

    public:
        /** \brief Constructor with a given range for each digit */
        MultiIndex(const array<int,dim>& limits)
            : limits_(limits) {
            std::fill(this->begin(), this->end(), 0);
        }

        /** \brief Increment the MultiIndex */
        MultiIndex& operator++() {

            for (int i=0; i<dim; i++) {

                // Augment digit
                (*this)[i]++;

                // If there is no carry-over we can stop here
                if ((*this)[i]<limits_[i])
                    break;

                (*this)[i] = 0;

            }
            return *this;
        }

        /** \brief Compute how many times you can call operator++ before getting to (0,...,0) again */
        int cycle() const {
            int result = 1;
            for (int i=0; i<dim; i++)
                result *= limits_[i];
            return result;
        }

    };

    /** \brief Insert a structured set of vertices into the factory */
    static void insertVertices(GridFactory<GridType>& factory,
                               const FieldVector<ctype,dim>& lowerLeft,
                               const FieldVector<ctype,dim>& upperRight,
                               const array<int,dim>& vertices) {

        // Pad lowerLeft and upperRight with zeros to dimworld dimensions
        FieldVector<ctype,dimworld> fullLowerLeft(0.);
        for (int i = 0; i < dim; ++i)
            fullLowerLeft[i] = lowerLeft[i];
        FieldVector<ctype,dimworld> fullUpperRight(0.);
        for (int i = 0; i < dim; ++i)
            fullUpperRight[i] = upperRight[i];

        MultiIndex index(vertices);

        // Compute the total number of vertices to be created
        int numVertices = index.cycle();

        // Create vertices
        for (int i=0; i<numVertices; i++, ++index) {

            // scale the multiindex to obtain a world position
            FieldVector<double,dimworld> pos(0);
            for (int j=0; j<dim; j++)
                pos[j] = fullLowerLeft[j] + index[j] * (fullUpperRight[j]-fullLowerLeft[j])/(vertices[j]-1);

            factory.insertVertex(pos);

        }

    }

    // Compute the index offsets needed to move to the adjacent vertices
    // in the different coordinate directions
    static array<int, dim> computeUnitOffsets(const array<int,dim>& vertices) {
        array<int, dim> unitOffsets;
        if (dim>0)  // paranoia
            unitOffsets[0] = 1;

        for (int i=1; i<dim; i++)
            unitOffsets[i] = unitOffsets[i-1] * vertices[i-1];

        return unitOffsets;
    }

public:

    /** \brief Create a structured cube grid
        \param lowerLeft Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements Number of elements in each coordinate direction
    */
    static std::auto_ptr<GridType> createCubeGrid(const FieldVector<ctype,dim>& lowerLeft,
            const FieldVector<ctype,dim>& upperRight,
            const array<int,dim>& elements) {
        // The grid factory
        GridFactory<GridType> factory;

        if (MPIHelper::getCollectiveCommunication().rank() == 0) {
            // Insert uniformly spaced vertices
            array<int,dim> vertices = elements;
            for( int i = 0; i < vertices.size(); ++i )
                vertices[i]++;

            // Insert vertices for structured grid into the factory
            insertVertices(factory, lowerLeft, upperRight, vertices);

            // Compute the index offsets needed to move to the adjacent
            // vertices in the different coordinate directions
            array<int, dim> unitOffsets =
                computeUnitOffsets(vertices);

            // Compute an element template (the cube at (0,...,0).  All
            // other cubes are constructed by moving this template around
            int nCorners = 1<<dim;

            std::vector<int> cornersTemplate(nCorners,0);

            for (int i=0; i<nCorners; i++)
                for (int j=0; j<dim; j++)
                    if ( i & (1<<j) )
                        cornersTemplate[i] += unitOffsets[j];

            // Insert elements
            MultiIndex index(elements);

            // Compute the total number of elementss to be created
            int numElements = index.cycle();

            for (int i=0; i<numElements; i++, ++index) {

                // 'base' is the index of the lower left element corner
                int base = 0;
                for (int j=0; j<dim; j++)
                    base += index[j] * unitOffsets[j];

                // insert new element
                std::vector<int> corners = cornersTemplate;
                for (size_t j=0; j<corners.size(); j++)
                    corners[j] += base;

                factory.insertElement
                (GeometryType(GeometryType::cube, dim), corners);

            }

        } // if(rank == 0)

        // Create the grid and hand it to the calling method
        return std::auto_ptr<GridType>(factory.createGrid());

    }

    /** \brief Create a structured simplex grid

        This works in all dimensions.  The Coxeter-Freudenthal-Kuhn triangulation is
        used, which splits each cube into dim! simplices.  See Allgower and Georg,
        'Numerical Path Following' for a description.
    */
    static std::auto_ptr<GridType> createSimplexGrid(const FieldVector<ctype,dim>& lowerLeft,
            const FieldVector<ctype,dim>& upperRight,
            const array<int,dim>& elements) {
        // The grid factory
        GridFactory<GridType> factory;

        if(MPIHelper::getCollectiveCommunication().rank() == 0) {
            // Insert uniformly spaced vertices
            array<int,dim> vertices = elements;
            for (std::size_t i=0; i<vertices.size(); i++)
                vertices[i]++;

            insertVertices(factory, lowerLeft, upperRight, vertices);

            // Compute the index offsets needed to move to the adjacent
            // vertices in the different coordinate directions
            array<int, dim> unitOffsets =
                computeUnitOffsets(vertices);

            // Insert the elements
            std::vector<int> corners(dim+1);

            // Loop over all "cubes", and split up each cube into dim!
            // (factorial) simplices
            MultiIndex elementsIndex(elements);
            int cycle = elementsIndex.cycle();

            for (int i=0; i<cycle; ++elementsIndex, i++) {

                // 'base' is the index of the lower left element corner
                int base = 0;
                for (int j=0; j<dim; j++)
                    base += elementsIndex[j] * unitOffsets[j];

                // each permutation of the unit vectors gives a simplex.
                std::vector<unsigned int> permutation(dim);
                for (int j=0; j<dim; j++)
                    permutation[j] = j;

                do {

                    // Make a simplex
                    std::vector<unsigned int> corners(dim+1);
                    corners[0] = base;

                    for (int j=0; j<dim; j++)
                        corners[j+1] =
                            corners[j] + unitOffsets[permutation[j]];

                    factory.insertElement
                    (GeometryType(GeometryType::simplex, dim),
                     corners);

                } while (std::next_permutation(permutation.begin(),
                                               permutation.end()));

            }

        } // if(rank == 0)

        // Create the grid and hand it to the calling method
        return std::auto_ptr<GridType>(factory.createGrid());
    }

};

/** \brief Specialization of the StructuredGridFactory for YaspGrid

    This allows a YaspGrid to be constructed using the
    StructuredGridFactory just like the unstructured Grids.  There are two
    limitations:
    \li YaspGrid does not support simplices
    \li YaspGrid only support grids which have their lower left corder at
        the origin.
*/
template<int dim>
class BemppStructuredGridFactory<YaspGrid<dim> >
{
    typedef YaspGrid<dim> GridType;
    typedef typename GridType::ctype ctype;
    static const int dimworld = GridType::dimensionworld;

public:
    /** \brief Create a structured cube grid

        \param lowerLeft  Lower left corner of the grid
        \param upperRight Upper right corner of the grid
        \param elements   Number of elements in each coordinate direction

        \note YaspGrid only supports lowerLeft at the origin.  This
              function throws a GridError if this requirement is not met.
    */
    static std::auto_ptr<GridType>
    createCubeGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                   const FieldVector<ctype,dimworld>& upperRight,
                   const array<int,dim>& elements) {
        for(int d = 0; d < dimworld; ++d)
            if(std::abs(lowerLeft[d]) > std::abs(upperRight[d])*1e-10)
                DUNE_THROW(GridError, className<BemppStructuredGridFactory>()
                           << "::createCubeGrid(): The lower coordinates "
                           "must be at the origin for YaspGrid.");

        FieldVector<int, dim> elements_;
        std::copy(elements.begin(), elements.end(), elements_.begin());

        return std::auto_ptr<GridType>
               (new GridType(upperRight, elements_,
                             FieldVector<bool,dim>(false), 0));
    }

    /** \brief Create a structured simplex grid

        \note Simplices are not supported in YaspGrid, so this functions
              unconditionally throws a GridError.
    */
    static std::auto_ptr<GridType>
    createSimplexGrid(const FieldVector<ctype,dimworld>& lowerLeft,
                      const FieldVector<ctype,dimworld>& upperRight,
                      const array<int,dim>& elements) {
        DUNE_THROW(GridError, className<BemppStructuredGridFactory>()
                   << "::createSimplexGrid(): Simplices are not supported "
                   "by YaspGrid.");
    }

};

}  // namespace Dune

#endif
