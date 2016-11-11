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

#include "cuda_grid.hpp"

#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <thrust/gather.h>
#include <thrust/transform.h>
#include <thrust/tabulate.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>

namespace Bempp {

  template <typename CoordinateType>
  struct calculateElementNormalAndIntegrationElementFunctor {

    __host__ __device__
    thrust::tuple<
        CoordinateType, CoordinateType, CoordinateType, CoordinateType>
    operator()(
        const thrust::tuple<CoordinateType, CoordinateType, CoordinateType,
                            CoordinateType, CoordinateType, CoordinateType,
                            CoordinateType, CoordinateType, CoordinateType>&
        elementCornerCoo) const {

      const CoordinateType vtx0x = thrust::get<0>(elementCornerCoo);
      const CoordinateType vtx0y = thrust::get<1>(elementCornerCoo);
      const CoordinateType vtx0z = thrust::get<2>(elementCornerCoo);

      const CoordinateType vtx1x = thrust::get<3>(elementCornerCoo);
      const CoordinateType vtx1y = thrust::get<4>(elementCornerCoo);
      const CoordinateType vtx1z = thrust::get<5>(elementCornerCoo);

      const CoordinateType vtx2x = thrust::get<6>(elementCornerCoo);
      const CoordinateType vtx2y = thrust::get<7>(elementCornerCoo);
      const CoordinateType vtx2z = thrust::get<8>(elementCornerCoo);

      CoordinateType nx = (vtx1y - vtx0y) * (vtx2z - vtx0z)
                - (vtx1z - vtx0z) * (vtx2y - vtx0y);

      CoordinateType ny = (vtx1z - vtx0z) * (vtx2x - vtx0x)
                - (vtx1x - vtx0x) * (vtx2z - vtx0z);

      CoordinateType nz = (vtx1x - vtx0x) * (vtx2y - vtx0y)
                - (vtx1y - vtx0y) * (vtx2x - vtx0x);

      const CoordinateType integrationElement = std::sqrt(nx*nx + ny*ny + nz*nz);

      nx /= integrationElement;
      ny /= integrationElement;
      nz /= integrationElement;

      return thrust::make_tuple(nx, ny, nz, integrationElement);
    }
  };

  template <typename CoordinateType>
  struct local2globalFunctor {

    unsigned int elemCount;

    thrust::device_ptr<const CoordinateType> vtx0x;
    thrust::device_ptr<const CoordinateType> vtx0y;
    thrust::device_ptr<const CoordinateType> vtx0z;

    thrust::device_ptr<const CoordinateType> vtx1x;
    thrust::device_ptr<const CoordinateType> vtx1y;
    thrust::device_ptr<const CoordinateType> vtx1z;

    thrust::device_ptr<const CoordinateType> vtx2x;
    thrust::device_ptr<const CoordinateType> vtx2y;
    thrust::device_ptr<const CoordinateType> vtx2z;

    thrust::device_ptr<const CoordinateType> fun0;
    thrust::device_ptr<const CoordinateType> fun1;
    thrust::device_ptr<const CoordinateType> fun2;

    local2globalFunctor(
      const unsigned int _elemCount,
      const thrust::device_ptr<const CoordinateType> _vtx0x,
      const thrust::device_ptr<const CoordinateType> _vtx0y,
      const thrust::device_ptr<const CoordinateType> _vtx0z,
      const thrust::device_ptr<const CoordinateType> _vtx1x,
      const thrust::device_ptr<const CoordinateType> _vtx1y,
      const thrust::device_ptr<const CoordinateType> _vtx1z,
      const thrust::device_ptr<const CoordinateType> _vtx2x,
      const thrust::device_ptr<const CoordinateType> _vtx2y,
      const thrust::device_ptr<const CoordinateType> _vtx2z,
      const thrust::device_ptr<const CoordinateType> _fun0,
      const thrust::device_ptr<const CoordinateType> _fun1,
      const thrust::device_ptr<const CoordinateType> _fun2)
      : elemCount(_elemCount),
        vtx0x(_vtx0x), vtx0y(_vtx0y), vtx0z(_vtx0z),
        vtx1x(_vtx1x), vtx1y(_vtx1y), vtx1z(_vtx1z),
        vtx2x(_vtx2x), vtx2y(_vtx2y), vtx2z(_vtx2z),
        fun0(_fun0), fun1(_fun1), fun2(_fun2) {}

    __host__ __device__
    thrust::tuple<CoordinateType, CoordinateType, CoordinateType> operator()(
        const unsigned int i) const {

      const unsigned int localPointIdx = i / elemCount;
      const unsigned int elementIdx = i % elemCount;

      const CoordinateType elVtx0x = vtx0x[elementIdx];
      const CoordinateType elVtx0y = vtx0y[elementIdx];
      const CoordinateType elVtx0z = vtx0z[elementIdx];

      const CoordinateType elVtx1x = vtx1x[elementIdx];
      const CoordinateType elVtx1y = vtx1y[elementIdx];
      const CoordinateType elVtx1z = vtx1z[elementIdx];

      const CoordinateType elVtx2x = vtx2x[elementIdx];
      const CoordinateType elVtx2y = vtx2y[elementIdx];
      const CoordinateType elVtx2z = vtx2z[elementIdx];

      const CoordinateType ptFun0 = fun0[localPointIdx];
      const CoordinateType ptFun1 = fun1[localPointIdx];
      const CoordinateType ptFun2 = fun2[localPointIdx];

      const CoordinateType xGlobal = ptFun0 * elVtx0x
                                   + ptFun1 * elVtx1x
                                   + ptFun2 * elVtx2x;
      const CoordinateType yGlobal = ptFun0 * elVtx0y
                                   + ptFun1 * elVtx1y
                                   + ptFun2 * elVtx2y;
      const CoordinateType zGlobal = ptFun0 * elVtx0z
                                   + ptFun1 * elVtx1z
                                   + ptFun2 * elVtx2z;

      return thrust::make_tuple(xGlobal, yGlobal, zGlobal);
    }
  };

  template <typename CoordinateType>
  CudaGrid<CoordinateType>::CudaGrid() {

    // Initialise member variables
    m_dim = 0;
    m_IdxCount = 0;
    m_VtxCount = 0;
    m_ElemCount = 0;

    m_vertices.clear();
    m_elementCorners.clear();

    m_vtx0x.clear();
    m_vtx0y.clear();
    m_vtx0z.clear();

    m_vtx1x.clear();
    m_vtx1y.clear();
    m_vtx1z.clear();

    m_vtx2x.clear();
    m_vtx2y.clear();
    m_vtx2z.clear();

    m_setupDone = false;
  }

  template <typename CoordinateType>
  CudaGrid<CoordinateType>::~CudaGrid() {

  }

  template <typename CoordinateType>
  void CudaGrid<CoordinateType>::pushGeometry(
      const Matrix<CoordinateType> &vertices,
      const Matrix<int> &elementCorners) {

    // Determine mesh parameters
    m_dim = vertices.cols();
    m_VtxCount = vertices.rows();
    m_IdxCount = elementCorners.cols();
    m_ElemCount = elementCorners.rows();

    if (m_dim != 3 || m_IdxCount != 3)
      throw std::runtime_error("CudaGrid::pushGeometry(): "
                               "only valid for triangular meshes in three dimensions.");

    // Allocate device memory
    m_vertices.resize(m_dim * m_VtxCount);
    m_elementCorners.resize(m_IdxCount * m_ElemCount);

    // Copy data to device
    m_vertices.assign(vertices.data(), vertices.data()+m_dim*m_VtxCount);
    m_elementCorners.assign(elementCorners.data(), elementCorners.data()+m_IdxCount*m_ElemCount);

    m_setupDone = false;

//    std::cout << "m_vertices = " << std::endl;
//    for (int i = 0; i < m_nVtx; ++i) {
//      for (int j = 0; j < m_dim; ++j) {
//        std::cout << m_vertices[j * m_nVtx + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    std::cout << "m_elementCorners = " << std::endl;
//    for (int i = 0; i < m_nEls; ++i) {
//      for (int j = 0; j < m_nIdx; ++j) {
//        std::cout << m_elementCorners[j * m_nEls + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
  }

  template <typename CoordinateType>
  void CudaGrid<CoordinateType>::setupGeometry() {

    if (m_setupDone == false) {

      // Allocate memory
      m_vtx0x.resize(m_ElemCount);
      m_vtx0y.resize(m_ElemCount);
      m_vtx0z.resize(m_ElemCount);

      m_vtx1x.resize(m_ElemCount);
      m_vtx1y.resize(m_ElemCount);
      m_vtx1z.resize(m_ElemCount);

      m_vtx2x.resize(m_ElemCount);
      m_vtx2y.resize(m_ElemCount);
      m_vtx2z.resize(m_ElemCount);

      // Measure time of the GPU execution (CUDA event based)
      cudaEvent_t start, stop;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);

      // Vertex 0
      thrust::gather(m_elementCorners.begin(),
                     m_elementCorners.begin()+m_ElemCount,
                     m_vertices.begin(),
                     m_vtx0x.begin());
      thrust::gather(m_elementCorners.begin(),
                     m_elementCorners.begin()+m_ElemCount,
                     m_vertices.begin()+m_VtxCount,
                     m_vtx0y.begin());
      thrust::gather(m_elementCorners.begin(),
                     m_elementCorners.begin()+m_ElemCount,
                     m_vertices.begin()+2*m_VtxCount,
                     m_vtx0z.begin());

      // Vertex 1
      thrust::gather(m_elementCorners.begin()+m_ElemCount,
                     m_elementCorners.begin()+2*m_ElemCount,
                     m_vertices.begin(),
                     m_vtx1x.begin());
      thrust::gather(m_elementCorners.begin()+m_ElemCount,
                     m_elementCorners.begin()+2*m_ElemCount,
                     m_vertices.begin()+m_VtxCount,
                     m_vtx1y.begin());
      thrust::gather(m_elementCorners.begin()+m_ElemCount,
                     m_elementCorners.begin()+2*m_ElemCount,
                     m_vertices.begin()+2*m_VtxCount,
                     m_vtx1z.begin());

      // Vertex 2
      thrust::gather(m_elementCorners.begin()+2*m_ElemCount,
                     m_elementCorners.end(),
                     m_vertices.begin(),
                     m_vtx2x.begin());
      thrust::gather(m_elementCorners.begin()+2*m_ElemCount,
                     m_elementCorners.end(),
                     m_vertices.begin()+m_VtxCount,
                     m_vtx2y.begin());
      thrust::gather(m_elementCorners.begin()+2*m_ElemCount,
                     m_elementCorners.end(),
                     m_vertices.begin()+2*m_VtxCount,
                     m_vtx2z.begin());

      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      float elapsedTimeGather;
      cudaEventElapsedTime(&elapsedTimeGather , start, stop);
      std::cout << "Time for gathering element corner coordinates is "
        << elapsedTimeGather << " ms" << std::endl;

      m_setupDone = true;

//        std::cout << "m_vtx0x = " << std::endl;
//        for (int i = 0; i < m_ElemCount; ++i) {
//          std::cout << m_vtx0x[i] << " " << std::flush;
//        }
//        std::cout << std::endl;
    }
  }

  template <typename CoordinateType>
  void CudaGrid<CoordinateType>::getRawGeometryData(
      unsigned int &vtxCount,
      unsigned int &elemCount,
      thrust::device_ptr<const CoordinateType> &vertices,
      thrust::device_ptr<const int> &elementCorners) {

    vtxCount = m_VtxCount;
    elemCount = m_ElemCount;
    vertices = m_vertices.data();
    elementCorners = m_elementCorners.data();
  }

  template <typename CoordinateType>
  void CudaGrid<CoordinateType>::calculateNormalsAndIntegrationElements(
      thrust::device_vector<CoordinateType> &normals,
      thrust::device_vector<CoordinateType> &integrationElements) const {

    if (m_setupDone == true) {

      // Allocate device memory
      normals.resize(m_dim * m_ElemCount);
      integrationElements.resize(m_ElemCount);

      // Measure time of the GPU execution (CUDA event based)
      cudaEvent_t start, stop;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);

      // Calculate element normals vectors and integration elements
      thrust::transform(
        thrust::make_zip_iterator(
          thrust::make_tuple(m_vtx0x.begin(), m_vtx0y.begin(), m_vtx0z.begin(),
                             m_vtx1x.begin(), m_vtx1y.begin(), m_vtx1z.begin(),
                             m_vtx2x.begin(), m_vtx2y.begin(), m_vtx2z.begin())),
        thrust::make_zip_iterator(
          thrust::make_tuple(m_vtx0x.end(), m_vtx0y.end(), m_vtx0z.end(),
                             m_vtx1x.end(), m_vtx1y.end(), m_vtx1z.end(),
                             m_vtx2x.end(), m_vtx2y.end(), m_vtx2z.end())),
        thrust::make_zip_iterator(
          thrust::make_tuple(normals.begin(), normals.begin()+m_ElemCount,
                             normals.begin()+2*m_ElemCount, integrationElements.begin())),
        calculateElementNormalAndIntegrationElementFunctor<CoordinateType>());

      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      float elapsedTimeNormals;
      cudaEventElapsedTime(&elapsedTimeNormals , start, stop);
      std::cout << "Time for calculating normals and integration elements is "
        << elapsedTimeNormals << " ms" << std::endl;

//    std::cout << "normals = " << std::endl;
//    for (int i = 0; i < m_ElemCount; ++i) {
//      for (int j = 0; j < m_dim; ++j) {
//        std::cout << normals[j * m_ElemCount + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//
//    std::cout << "integrationElements = " << std::endl;
//    for (int i = 0; i < m_ElemCount; ++i) {
//      std::cout << integrationElements[i] << " " << std::flush;
//    }
//    std::cout << std::endl;

    } else {
      throw std::runtime_error("CudaGrid::calculateNormalsAndIntegrationElements(): "
                               "setup required");
    }
  }

  template <typename CoordinateType>
  void CudaGrid<CoordinateType>::local2global(
      const Matrix<CoordinateType> &localPoints,
      thrust::device_vector<CoordinateType> &globalPoints) const {

    if (m_setupDone == true) {

      const unsigned int localPointDim = localPoints.rows();
      const unsigned int localPointCount = localPoints.cols();

      if (localPointDim != 2)
        throw std::runtime_error("CudaGrid::local2global(): "
                                 "only valid for two-dimensional local points");

      // Allocate device memory
      globalPoints.resize(m_dim * m_ElemCount * localPointCount);
      // [xPt0el0 xPt0el1 ... xPt0elM | xPt1el0 ... xPt1elM | ... | ... xPtNelM |
      //  yPt0el0 yPt0el1 ... yPt0elM | yPt1el0 ... yPt1elM | ... | ... yPtNelM |
      //  zPt0el0 zPt0el1 ... zPt0elM | zPt1el0 ... zPt1elM | ... | ... zPtNelM ]

      thrust::host_vector<CoordinateType> h_geomShapeFun0(localPointCount);
      thrust::host_vector<CoordinateType> h_geomShapeFun1(localPointCount);
      thrust::host_vector<CoordinateType> h_geomShapeFun2(localPointCount);
      // [pt0 pt1 ... ptN]

      // Evaluate geometrical shape function values
      for (int localPoint = 0; localPoint < localPointCount; ++localPoint) {
        const CoordinateType r = localPoints(0,localPoint);
        const CoordinateType s = localPoints(1,localPoint);
        h_geomShapeFun0[localPoint] = 1.0 - r - s;
        h_geomShapeFun1[localPoint] = r;
        h_geomShapeFun2[localPoint] = s;
      }

      // Copy data to device
      thrust::device_vector<CoordinateType> d_geomShapeFun0 = h_geomShapeFun0;
      thrust::device_vector<CoordinateType> d_geomShapeFun1 = h_geomShapeFun1;
      thrust::device_vector<CoordinateType> d_geomShapeFun2 = h_geomShapeFun2;

//    std::cout << "d_geomShapeFun = " << std::endl;
//    for (int i = 0; i < localPointCount; ++i) {
//      std::cout << d_geomShapeFun0[i] << " "
//                << d_geomShapeFun1[i] << " "
//                << d_geomShapeFun2[i] << std::endl;
//    }
//    std::cout << std::endl;

      // Measure time of the GPU execution (CUDA event based)
      cudaEvent_t start, stop;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start, 0);

      thrust::tabulate(
        thrust::make_zip_iterator(
          thrust::make_tuple(globalPoints.begin(),
                             globalPoints.begin()+m_ElemCount*localPointCount,
                             globalPoints.begin()+2*m_ElemCount*localPointCount)),
        thrust::make_zip_iterator(
          thrust::make_tuple(globalPoints.begin()+m_ElemCount*localPointCount,
                             globalPoints.begin()+2*m_ElemCount*localPointCount,
                             globalPoints.end())),
        local2globalFunctor<CoordinateType>(m_ElemCount,
                            m_vtx0x.data(), m_vtx0y.data(), m_vtx0z.data(),
                            m_vtx1x.data(), m_vtx1y.data(), m_vtx1z.data(),
                            m_vtx2x.data(), m_vtx2y.data(), m_vtx2z.data(),
                            d_geomShapeFun0.data(),
                            d_geomShapeFun1.data(),
                            d_geomShapeFun2.data()));

      cudaEventRecord(stop, 0);
      cudaEventSynchronize(stop);
      float elapsedTimeMapping;
      cudaEventElapsedTime(&elapsedTimeMapping , start, stop);
      std::cout << "Time for mapping local to global coordinates is "
        << elapsedTimeMapping << " ms" << std::endl;

//    std::cout << "globalPoints = " << std::endl;
//    for (int i = 0; i < activeElemCount; ++i) {
//      for (int j = 0; j < localPointCount; ++j) {
//        for (int k = 0; k < m_dim; ++k) {
//          std::cout << globalPoints[k * localPointCount * activeElemCount
//                                  + j * activeElemCount
//                                  + i] << " " << std::flush;
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//    }

    } else {
      throw std::runtime_error("CudaGrid::local2globalImpl(): "
                               "setup required");
    }
  }

  // TODO: Instantiate class templared on coordinate type?
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_DP_REAL(CudaGrid);
  FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_SP_REAL(CudaGrid);

} // namespace Bempp
