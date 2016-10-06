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

#include "cuda_handler.hpp"

#ifdef WITH_CUDA

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace Fiber {

CUDAHandler::CUDAHandler(const bool useCUDA, const CUDAOptions &options) {

  this->useCUDA = useCUDA;

  deviceUsed = options.deviceUsed;
}

CUDAHandler::~CUDAHandler() {

}

template <typename CoordinateType, typename IndexType>
void CUDAHandler::pushGeometry(const Matrix<CoordinateType> &vtx,
                               const Matrix<IndexType> &idx) const {
  // Allocate memory on device
  meshGeom.size.dim = vtx.rows();
  meshGeom.size.nVtx = vtx.cols();
  size_t vtxSize = meshGeom.size.dim * meshGeom.size.nVtx;
  thrust::device_vector<CoordinateType> cu_vtx(vtxSize);

  meshGeom.size.nEls = idx.cols();
  meshGeom.size.nIdx = idx.rows();
  size_t idxSize = meshGeom.size.nEls * meshGeom.size.nIdx;
  thrust::device_vector<IndexType> cu_idx(idxSize);

  // TODO: Copy the geometry data to device memory

}

} // namespace Fiber

#endif // WITH_CUDA
