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
   * \param[in] vertices mesh vertices
   * \param[in] elementCorners element corner indices
   */
  void pushGeometry(const Matrix<double> &vertices,
                    const Matrix<int> &elementCorners);

  /**
   * \brief Setup geometry data for specified elements on the device i.e. gather
   * element corner coordinates, calculate normal vectors and integration elements
   * \param[in] elementIndices indices of elements to set up
   */
  void setupElements(const std::vector<int> &elementIndices);

  /**
   * \brief Clean up element data on the device
   */
  void freeElementData();

  /**
   * \brief Convert local (logical) to global (physical) coordinates on the device
   * \param[in] localPoints local coordinates
   * \param[out] globalPoints global coordinates
   */
  void local2global(const Matrix<double> &localPoints,
                    thrust::device_vector<double> &globalPoints) const;

private:
  /** \cond PRIVATE */

  /**
   * \brief Setup geometry data for all elements on the device, i.e. gather
   * element corner coordinates, calculate normal vectors and integration elements
   * */
  void setupGeometry();

  /**
   * \brief Calculate element normal vectors and determinants of Jacobian
   * (factor appearing in the integral transformation formula) on the device
   */
  void calculateNormalsAndIntegrationElements();

  // Mesh parameters
  unsigned int m_dim;
  unsigned int m_IdxCount;
  unsigned int m_VtxCount;
  unsigned int m_ElemCount;

  thrust::device_vector<double> m_vertices;
  thrust::device_vector<int> m_elementCorners;

  // Element corner coordinates
  thrust::device_vector<double> m_vtx0x;
  thrust::device_vector<double> m_vtx0y;
  thrust::device_vector<double> m_vtx0z;

  thrust::device_vector<double> m_vtx1x;
  thrust::device_vector<double> m_vtx1y;
  thrust::device_vector<double> m_vtx1z;

  thrust::device_vector<double> m_vtx2x;
  thrust::device_vector<double> m_vtx2y;
  thrust::device_vector<double> m_vtx2z;

  thrust::device_vector<double> m_normals;
  thrust::device_vector<double> m_integrationElements;

  thrust::device_vector<int> m_activeElemIndices;

  /** \endcond */
};

} // namespace Bempp

#endif
