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

#ifndef bempp_grid_factory_hpp
#define bempp_grid_factory_hpp

#include "common.hpp"

#include <armadillo>
#include <memory>

namespace Bempp
{

// Forward declarations
class Grid;

/** \brief %Grid parameters.

  This structure is used to specify parameters of grid constructed by GridFactory.
  */
struct GridParameters {
    /** \brief %Grid topology */
    enum Topology {
        LINEAR, /**< \brief one-dimensional grid */
        TRIANGULAR, /**< \brief grid composed of triangular elements */
        QUADRILATERAL, /**< \brief grid composed of quadrilateral elements */
        HYBRID /**< \brief grid composed (potentially) both of triangular and quadrilateral elements */
    } topology;
};

/** \brief %Grid factory.

  This class provides static member functions to construct grids on the fly and to import grids from
  existing files.
  */
class GridFactory
{
public:
    /** \brief Construct a regular structured grid.

      \param params     Parameters of the grid to be constructed.
      \param lowerLeft  Coordinates of the lower left corner of the grid.
      \param upperRight Coordinates of the upper right corner of the grid.
      \param nElements  Number of grid subdivisions in each direction.

      This function constructs a regular structured grid. Its dimension, \p
      dimGrid, and the dimension of the surrounding space, \p dimWorld, are determined from the parameter \p params.
      The constructed grid covers the \p dimGrid-dimensional cube

      [lowerLeft(0) upperRight(0)] x [lowerLeft(1) upperRight(1)] x ... x [lowerLeft(dimGrid-1), upperRight(dimGrid-1)].

      The last \p dimWorld - \p dimGrid dimensions of all grid points are set to zero.

      Each side of the cube parallel to the <em>n</em>th coordinate axis is subdivided into nElements(n) segments.

      \note Currently only grids with trangular topology are supported.
    */
    static std::auto_ptr<Grid> createStructuredGrid(const GridParameters& params,
            const arma::Col<ctype>& lowerLeft,
            const arma::Col<ctype>& upperRight,
            const arma::Col<unsigned int>& nElements);

    /** \brief Import grid from a file in Gmsh format.

      \param params Parameters of the grid to be constructed.
      \param fileName Name of the Gmsh file.
      \param verbose  Output diagnostic information.
      \param insertBoundarySegments

      \bug Ask Dune developers about the significance of insertBoundarySegments.
      \see <a href>http://geuz.org/gmsh/</a> for information about the Gmsh file format.
      \see Dune::GmshReader documentation for information about the supported Gmsh features.
    */
    static std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
            const std::string& fileName,
            bool verbose=true, bool insertBoundarySegments=true);

    /** \brief Import grid from a file in Gmsh format.

      \param params Parameters of the grid to be constructed.
      \param fileName Name of the Gmsh file.
      \param boundaryId2PhysicalEntity
      \param elementIndex2PhysicalEntity
      \param verbose  Output diagnostic information.
      \param insertBoundarySegments

      \bug Ask Dune developers about the significance of the undocumented parameters.
      \see <a href>http://geuz.org/gmsh/</a> for information about the Gmsh file format.
      \see Dune::GmshReader documentation for information about the supported Gmsh features.
    */
    static std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
            const std::string& fileName,
            std::vector<int>& boundaryId2PhysicalEntity,
            std::vector<int>& elementIndex2PhysicalEntity,
            bool verbose=true, bool insertBoundarySegments=true);
};

} // namespace Bempp

#endif
