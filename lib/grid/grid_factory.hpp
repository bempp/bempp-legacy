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

#ifndef bempp_grid_factory_hpp
#define bempp_grid_factory_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "grid_parameters.hpp"
#include "../common/eigen_support.hpp"

#include <memory>

namespace Bempp {

/** \cond FORWARD_DECL */
class Grid;
/** \endcond */

/** \ingroup grid
    \brief %Grid factory.

  This class provides static member functions to construct grids on the fly and
  to import grids from
  existing files.
  */
class GridFactory {
public:
  /** \brief Construct a regular structured grid.

    \param[in] params     Parameters of the grid to be constructed.
    \param[in] lowerLeft  Coordinates of the lower left corner of the grid.
    \param[in] upperRight Coordinates of the upper right corner of the grid.
    \param[in] nElements  Number of grid subdivisions in each direction.

    This function constructs a regular structured grid. Its dimension, \p
    dimGrid, and the dimension of the surrounding space, \p dimWorld, are
    determined from the parameter \p params. The constructed grid covers the
    \p dimGrid-dimensional cube

    [lowerLeft(0) upperRight(0)] x [lowerLeft(1) upperRight(1)] x ... x
    [lowerLeft(dimGrid-1), upperRight(dimGrid-1)].

    The last \p dimWorld - \p dimGrid dimensions of all grid points are set
    to zero.

    Each side of the cube parallel to the <em>n</em>th coordinate axis is
    subdivided into nElements(n) segments.

    \note Currently only grids with triangular topology are supported.
  */
  static shared_ptr<Grid>
  createStructuredGrid(const GridParameters &params,
                       const Vector<double> &lowerLeft,
                       const Vector<double> &upperRight,
                       const Vector<unsigned int> &nElements);

  /** \brief Import grid from a file in Gmsh format.

    \param[in] params Parameters of the grid to be constructed.
    \param[in] fileName Name of the Gmsh file.
    \param[in] verbose  Output diagnostic information.
    \param[in] insertBoundarySegments

    \bug Ask Dune developers about the significance of insertBoundarySegments.
    \see <a href>http://geuz.org/gmsh/</a> for information about the Gmsh file
    format.
    \see Dune::GmshReader documentation for information about the supported Gmsh
    features.
  */
  static shared_ptr<Grid> importGmshGrid(const GridParameters &params,
                                         const std::string &fileName,
                                         bool verbose = false,
                                         bool insertBoundarySegments = false);

  /** \brief Import grid from a file in Gmsh format.

    \param[in] params Parameters of the grid to be constructed.
    \param[in] fileName Name of the Gmsh file.
    \param[in] boundaryId2PhysicalEntity
    \param[in] elementIndex2PhysicalEntity
    \param[in] verbose  Output diagnostic information.
    \param[in] insertBoundarySegments

    \bug Ask Dune developers about the significance of the undocumented
    parameters.
    \see <a href>http://geuz.org/gmsh/</a> for information about the Gmsh file
    format.
    \see Dune::GmshReader documentation for information about the supported Gmsh
    features.
  */
  static shared_ptr<Grid>
  importGmshGrid(const GridParameters &params, const std::string &fileName,
                 std::vector<int> &boundaryId2PhysicalEntity,
                 std::vector<int> &elementIndex2PhysicalEntity,
                 bool verbose = true, bool insertBoundarySegments = false);

  /** \brief Create a grid from connectivity arrays.
   *
   *  \param[in] params
   *    Parameters of the grid to be constructed.
   *  \param[in] vertices
   *    2D array whose (i, j)th element contains the ith
   *    component of the jth vertex.
   *  \param[in] elementCorners
   *    2D array whose (i, j)th element contains the index of the ith vertex
   *    of the jth element.
   *  \param[in] domainIndices
   *    (Optional) Vector whose ith element contains the domain index of
   *    the ith element. By default, this argument is set to an empty vector,
   *    in which case all elements are taken to belong to domain 0.
   *
   *  \note Currently only grids with triangular topology are supported.
   */
  static shared_ptr<Grid> createGridFromConnectivityArrays(
      const GridParameters &params, const Matrix<double> &vertices,
      const Matrix<int> &elementCorners,
      const std::vector<int> &domainIndices = std::vector<int>());
};

} // namespace Bempp

#endif
