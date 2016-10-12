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

#ifndef bempp_cuda_grid_hpp
#define bempp_cuda_grid_hpp

/** \file . */

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include <thrust/device_vector.h>

namespace Bempp {

class CudaGrid {
public:

  /** \brief Constructor */
  CudaGrid();

  /** \brief Destructor */
  virtual ~CudaGrid();

  /**
   * \brief Push the mesh geometry to device memory
   * \param vertices mesh vertices
   * \param elementCorners element corner indices
   */
  void pushGeometry(const Matrix<double> &vertices,
                    const Matrix<int> &elementCorners);
  /**
   * \brief Setup the geometry on the device, i.e. gather element corner
   * coordinates, calculate normal vectors and integration elements
   * */
  void setupGeometry();

  /**
   * \brief Calculate element normal vectors and determinants of Jacobian
   * (factor appearing in the integral transformation formula) on the device
   */
  void calculateNormalsAndIntegrationElements();

  /**
   * \brief Convert local (logical) to global (physical) coordinates for all
   * elements on the device
   * \param localPoints local coordinates
   * \param globalPoints global coordinates
   */
  void local2global(const thrust::host_vector<double> &localPoints,
                    thrust::device_vector<double> &globalPoints);

private:
  /** \cond PRIVATE */

  // Mesh parameters
  unsigned int dim;
  unsigned int nIdx;
  unsigned int nVtx;
  unsigned int nEls;

  thrust::device_vector<double> d_vertices;
  // [x0 x1 x2 ... xN | y0 y1 ... yN | z0 z1 ... zN]

  thrust::device_vector<int> d_elementCorners;
  // [vtx0el0 vtx0el1 ... vtx0elN | vtx1el0 ... vtx1elN | vtx2el0 ... vtx2elN]

  // Element corner coordinates
  thrust::device_vector<double> d_vtx0x;
  thrust::device_vector<double> d_vtx0y;
  thrust::device_vector<double> d_vtx0z;

  thrust::device_vector<double> d_vtx1x;
  thrust::device_vector<double> d_vtx1y;
  thrust::device_vector<double> d_vtx1z;

  thrust::device_vector<double> d_vtx2x;
  thrust::device_vector<double> d_vtx2y;
  thrust::device_vector<double> d_vtx2z;
  // [el0 el1 el2 ... elN]

  thrust::device_vector<double> d_normals;
  // [x0 x1 x2 ... xN | y0 y1 ... yN | z0 z1 ... zN]

  thrust::device_vector<double> d_integrationElements;
  // [el0 el1 el2 ... elN]

  /** \endcond */
};

} // namespace Bempp

#endif
