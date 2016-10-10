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

private:
  /** \cond PRIVATE */
  // TODO: double, int -> CoordinateType, IndexType
  thrust::device_vector<double> d_vertices;
  thrust::device_vector<int> d_elementCorners;
  thrust::device_vector<double> d_normals;
  thrust::device_vector<double> d_integrationElements;
  /** \endcond */
};

} // namespace Bempp

#endif
