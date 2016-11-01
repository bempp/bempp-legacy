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
   * \brief Setup geometry data for all elements on the device i.e. gather
   * element corner coordinates
   */
  void setupGeometry();

  /**
   * \brief Calculate element normal vectors and determinants of Jacobian
   * (factor appearing in the integral transformation formula) on the device
   */
  void calculateNormalsAndIntegrationElements(
      thrust::device_vector<double> &normals,
      thrust::device_vector<double> &integrationElements) const;
  /** \overload */
  void calculateNormalsAndIntegrationElements(
      thrust::device_vector<float> &normals,
      thrust::device_vector<float> &integrationElements) const;

  void getRawGeometryData(unsigned int &vtxCount,
                          unsigned int &elemCount,
                          thrust::device_ptr<const double> &vertices,
                          thrust::device_ptr<const int> &elementCorners);

  /**
   * \brief Convert local (logical) to global (physical) coordinates on the device
   * \param[in] localPoints Matrix whose \f$i\f$th column contains the
      local coordinates of a point \f$x_i \in D\f$
   * \param[out] globalPoints Vector containing the global coordinates
   */
  void local2global(const Matrix<double> &localPoints,
                    thrust::device_vector<double> &globalPoints) const;
  /** \overload */
  void local2global(const Matrix<float> &localPoints,
                    thrust::device_vector<float> &globalPoints) const;

private:
  /** \cond PRIVATE */

  void calculateNormalsAndIntegrationElementsImpl(
      thrust::device_vector<double> &normals,
      thrust::device_vector<double> &integrationElements) const;

  void local2globalImpl(const Matrix<double> &localPoints,
                        thrust::device_vector<double> &globalPoints) const;

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

  bool m_setupDone;

  // Helper functions for implementation
  template <typename T1, typename T2>
  void convertMat(const Matrix<T1> &in, Matrix<T2> &out) const;

  /** \endcond */
};

} // namespace Bempp

#endif
