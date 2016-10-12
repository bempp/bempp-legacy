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

#include <thrust/gather.h>
#include <thrust/transform.h>
#include <thrust/tabulate.h>
#include <thrust/execution_policy.h>

namespace Bempp {

  struct calculateElementNormalAndIntegrationElementFunctor {

    __host__ __device__
    thrust::tuple<double, double, double, double> operator()(
      const thrust::tuple<double, double, double,
                          double, double, double,
                          double, double, double>& elementCornerCoo) const {

      double vtx0x = thrust::get<0>(elementCornerCoo);
      double vtx0y = thrust::get<1>(elementCornerCoo);
      double vtx0z = thrust::get<2>(elementCornerCoo);

      double vtx1x = thrust::get<3>(elementCornerCoo);
      double vtx1y = thrust::get<4>(elementCornerCoo);
      double vtx1z = thrust::get<5>(elementCornerCoo);

      double vtx2x = thrust::get<6>(elementCornerCoo);
      double vtx2y = thrust::get<7>(elementCornerCoo);
      double vtx2z = thrust::get<8>(elementCornerCoo);

      double nx = (vtx1y - vtx0y) * (vtx2z - vtx0z)
                - (vtx1z - vtx0z) * (vtx2y - vtx0y);

      double ny = (vtx1z - vtx0z) * (vtx2x - vtx0x)
                - (vtx1x - vtx0x) * (vtx2z - vtx0z);

      double nz = (vtx1x - vtx0x) * (vtx2y - vtx0y)
                - (vtx1y - vtx0y) * (vtx2x - vtx0x);

      double integrationElement = std::sqrt(nx*nx + ny*ny + nz*nz);

      nx /= integrationElement;
      ny /= integrationElement;
      nz /= integrationElement;

      return thrust::make_tuple(nx, ny, nz, integrationElement);
    }
  };

  struct geomShapeFunFunctor {

    __host__ __device__
    thrust::tuple<double, double, double> operator()(
      const thrust::tuple<double, double>& localPointCoo) const {

      double r = thrust::get<0>(localPointCoo);
      double s = thrust::get<1>(localPointCoo);

      double fun0 = 1.0 - r - s;
      double fun1 = r;
      double fun2 = s;

      return thrust::make_tuple(fun0, fun1, fun2);
    }
  };

  struct local2globalFunctor {

    unsigned int nEls;

    thrust::device_ptr<double> vtx0x;
    thrust::device_ptr<double> vtx0y;
    thrust::device_ptr<double> vtx0z;

    thrust::device_ptr<double> vtx1x;
    thrust::device_ptr<double> vtx1y;
    thrust::device_ptr<double> vtx1z;

    thrust::device_ptr<double> vtx2x;
    thrust::device_ptr<double> vtx2y;
    thrust::device_ptr<double> vtx2z;

    thrust::device_ptr<double> fun0;
    thrust::device_ptr<double> fun1;
    thrust::device_ptr<double> fun2;

    local2globalFunctor(
      const unsigned int _nEls,
      thrust::device_ptr<double> _vtx0x,
      thrust::device_ptr<double> _vtx0y,
      thrust::device_ptr<double> _vtx0z,
      thrust::device_ptr<double> _vtx1x,
      thrust::device_ptr<double> _vtx1y,
      thrust::device_ptr<double> _vtx1z,
      thrust::device_ptr<double> _vtx2x,
      thrust::device_ptr<double> _vtx2y,
      thrust::device_ptr<double> _vtx2z,
      thrust::device_ptr<double> _fun0,
      thrust::device_ptr<double> _fun1,
      thrust::device_ptr<double> _fun2)
      : nEls(_nEls),
        vtx0x(_vtx0x), vtx0y(_vtx0y), vtx0z(_vtx0z),
        vtx1x(_vtx1x), vtx1y(_vtx1y), vtx1z(_vtx1z),
        vtx2x(_vtx2x), vtx2y(_vtx2y), vtx2z(_vtx2z),
        fun0(_fun0), fun1(_fun1), fun2(_fun2) {}

    __host__ __device__
    thrust::tuple<double, double, double> operator()(
        const unsigned int i) const {

      // Which memory mapping to go for?
//      unsigned int localPointIdx = i % nLocalPoints;
//      unsigned int elementIdx = i / nLocalPoints;
      unsigned int localPointIdx = i / nEls;
      unsigned int elementIdx = i % nEls;

      double elVtx0x = vtx0x[elementIdx];
      double elVtx0y = vtx0y[elementIdx];
      double elVtx0z = vtx0z[elementIdx];

      double elVtx1x = vtx1x[elementIdx];
      double elVtx1y = vtx1y[elementIdx];
      double elVtx1z = vtx1z[elementIdx];

      double elVtx2x = vtx2x[elementIdx];
      double elVtx2y = vtx2y[elementIdx];
      double elVtx2z = vtx2z[elementIdx];

      double ptFun0 = fun0[localPointIdx];
      double ptFun1 = fun1[localPointIdx];
      double ptFun2 = fun2[localPointIdx];

      double xGlobal = ptFun0 * elVtx0x
                     + ptFun1 * elVtx1x
                     + ptFun2 * elVtx2x;
      double yGlobal = ptFun0 * elVtx0y
                     + ptFun1 * elVtx1y
                     + ptFun2 * elVtx2y;
      double zGlobal = ptFun0 * elVtx0z
                     + ptFun1 * elVtx1z
                     + ptFun2 * elVtx2z;

      return thrust::make_tuple(xGlobal, yGlobal, zGlobal);
    }
  };

  CudaGrid::CudaGrid() {

    // Initialise member variables
    dim = 0;
    nIdx = 0;
    nVtx = 0;
    nEls = 0;

    d_vertices.clear();
    d_elementCorners.clear();

    d_vtx0x.clear();
    d_vtx0y.clear();
    d_vtx0z.clear();

    d_vtx1x.clear();
    d_vtx1y.clear();
    d_vtx1z.clear();

    d_vtx2x.clear();
    d_vtx2y.clear();
    d_vtx2z.clear();

    d_normals.clear();
    d_integrationElements.clear();
  }

  CudaGrid::~CudaGrid() {

  }

  void CudaGrid::pushGeometry(const Matrix<double> &vertices,
                              const Matrix<int> &elementCorners) {

    // Determine mesh parameters
    dim = vertices.cols();
    nVtx = vertices.rows();
    nIdx = elementCorners.cols();
    nEls = elementCorners.rows();

    if (dim != 3 || nIdx != 3)
      throw std::runtime_error("CudaGrid::pushGeometry(): "
                               "only valid for triangular meshes in three dimensions.");

    // Allocate device memory
    d_vertices.resize(dim * nVtx);
    d_elementCorners.resize(nIdx * nEls);

    // Copy data to device
    d_vertices.assign(vertices.data(), vertices.data()+dim*nVtx);
    d_elementCorners.assign(elementCorners.data(), elementCorners.data()+nIdx*nEls);

//    std::cout << "d_vertices = " << std::endl;
//    for (int i = 0; i < nVtx; ++i) {
//      for (int j = 0; j < dim; ++j) {
//        std::cout << d_vertices[j * nVtx + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    std::cout << "d_elementCorners = " << std::endl;
//    for (int i = 0; i < nEls; ++i) {
//      for (int j = 0; j < nIdx; ++j) {
//        std::cout << d_elementCorners[j * nEls + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;

    // TEST Setup geometry
    setupGeometry();
  }

  void CudaGrid::setupGeometry() {

    // Gather element corner coordinates
    d_vtx0x.resize(nEls);
    d_vtx0y.resize(nEls);
    d_vtx0z.resize(nEls);

    d_vtx1x.resize(nEls);
    d_vtx1y.resize(nEls);
    d_vtx1z.resize(nEls);

    d_vtx2x.resize(nEls);
    d_vtx2y.resize(nEls);
    d_vtx2z.resize(nEls);

    // Measure time of the GPU execution (CUDA event based)
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    // Vertex 0
    thrust::gather(d_elementCorners.begin(), d_elementCorners.begin()+nEls,
                   d_vertices.begin(),
                   d_vtx0x.begin());
    thrust::gather(d_elementCorners.begin(), d_elementCorners.begin()+nEls,
                   d_vertices.begin()+nVtx,
                   d_vtx0y.begin());
    thrust::gather(d_elementCorners.begin(), d_elementCorners.begin()+nEls,
                   d_vertices.begin()+2*nVtx,
                   d_vtx0z.begin());

    // Vertex 1
    thrust::gather(d_elementCorners.begin()+nEls, d_elementCorners.begin()+2*nEls,
                   d_vertices.begin(),
                   d_vtx1x.begin());
    thrust::gather(d_elementCorners.begin()+nEls, d_elementCorners.begin()+2*nEls,
                   d_vertices.begin()+nVtx,
                   d_vtx1y.begin());
    thrust::gather(d_elementCorners.begin()+nEls, d_elementCorners.begin()+2*nEls,
                   d_vertices.begin()+2*nVtx,
                   d_vtx1z.begin());

    // Vertex 2
    thrust::gather(d_elementCorners.begin()+2*nEls, d_elementCorners.end(),
                   d_vertices.begin(),
                   d_vtx2x.begin());
    thrust::gather(d_elementCorners.begin()+2*nEls, d_elementCorners.end(),
                   d_vertices.begin()+nVtx,
                   d_vtx2y.begin());
    thrust::gather(d_elementCorners.begin()+2*nEls, d_elementCorners.end(),
                   d_vertices.begin()+2*nVtx,
                   d_vtx2z.begin());

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeGather;
    cudaEventElapsedTime(&elapsedTimeGather , start, stop);
    std::cout << "Time for gathering element corner coordinates is "
      << elapsedTimeGather << " ms" << std::endl;

//    std::cout << "d_vtx0x = " << std::endl;
//    for (int i = 0; i < nEls; ++i) {
//      std::cout << d_vtx0x[i] << std::endl;
//    }
//    std::cout << std::endl;

    calculateNormalsAndIntegrationElements();

    // TEST Convert local to global coordinates for all elements
    thrust::host_vector<double> localPoints(2);
    localPoints[0] = 0.1;
    localPoints[1] = 0.1;

    thrust::device_vector<double> globalPoints;
    local2global(localPoints, globalPoints);
  }

  void CudaGrid::calculateNormalsAndIntegrationElements() {

    // Allocate device memory
    d_normals.resize(dim * nEls);
    d_integrationElements.resize(nEls);

    // Measure time of the GPU execution (CUDA event based)
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    // Calculate element normals vectors and integration elements
    thrust::transform(
      thrust::make_zip_iterator(
        thrust::make_tuple(d_vtx0x.begin(), d_vtx0y.begin(), d_vtx0z.begin(),
                           d_vtx1x.begin(), d_vtx1y.begin(), d_vtx1z.begin(),
                           d_vtx2x.begin(), d_vtx2y.begin(), d_vtx2z.begin())),
      thrust::make_zip_iterator(
        thrust::make_tuple(d_vtx0x.end(), d_vtx0y.end(), d_vtx0z.end(),
                           d_vtx1x.end(), d_vtx1y.end(), d_vtx1z.end(),
                           d_vtx2x.end(), d_vtx2y.end(), d_vtx2z.end())),
      thrust::make_zip_iterator(
        thrust::make_tuple(d_normals.begin(), d_normals.begin()+nEls,
                           d_normals.begin()+2*nEls, d_integrationElements.begin())),
      calculateElementNormalAndIntegrationElementFunctor());

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeNormals;
    cudaEventElapsedTime(&elapsedTimeNormals , start, stop);
    std::cout << "Time for calculating normals and integration elements is "
      << elapsedTimeNormals << " ms" << std::endl;

//    std::cout << "d_normals = " << std::endl;
//    for (int i = 0; i < nEls; ++i) {
//      for (int j = 0; j < dim; ++j) {
//        std::cout << d_normals[j * nEls + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    std::cout << "d_integrationElements = " << std::endl;
//    for (int i = 0; i < nEls; ++i) {
//      std::cout << d_integrationElements[i] << std::endl;
//    }
//    std::cout << std::endl;
  }

  void CudaGrid::local2global(const thrust::host_vector<double> &localPoints,
                              thrust::device_vector<double> &globalPoints) {

    // localPoints = [r0 r1 ... rN | s0 s1 ... sN]

    const unsigned int localPointDim = 2;
    const unsigned int nLocalPoints = localPoints.size() / localPointDim;

    if (localPoints.size() % localPointDim != 0)
      throw std::runtime_error("CudaGrid::local2global(): "
                               "only valid for two-dimensional local points");

    // Allocate device memory
    globalPoints.resize(dim * nEls * nLocalPoints);
    // [xPt0el0 xPt0el1 ... xPt0elM | xPt1el0 ... xPt1elM | ... | ... xPtNelM |
    //  yPt0el0 yPt0el1 ... yPt0elM | yPt1el0 ... yPt1elM | ... | ... yPtNelM |
    //  zPt0el0 zPt0el1 ... zPt0elM | zPt1el0 ... zPt1elM | ... | ... zPtNelM ]

    thrust::host_vector<double> h_geomShapeFun0(nLocalPoints);
    thrust::host_vector<double> h_geomShapeFun1(nLocalPoints);
    thrust::host_vector<double> h_geomShapeFun2(nLocalPoints);
    // [pt0 pt1 ... ptN]

    // Evaluate geometrical shape function values on the host
    thrust::transform(thrust::host,
      thrust::make_zip_iterator(
        thrust::make_tuple(localPoints.begin(),
                           localPoints.begin()+nLocalPoints)),
      thrust::make_zip_iterator(
        thrust::make_tuple(localPoints.begin()+nLocalPoints,
                           localPoints.end())),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_geomShapeFun0.begin(),
                           h_geomShapeFun1.begin(),
                           h_geomShapeFun2.begin())),
      geomShapeFunFunctor());

    // Copy data to device
    thrust::device_vector<double> geomShapeFun0 = h_geomShapeFun0;
    thrust::device_vector<double> geomShapeFun1 = h_geomShapeFun1;
    thrust::device_vector<double> geomShapeFun2 = h_geomShapeFun2;

//    std::cout << "geomShapeFun = " << std::endl;
//    for (int i = 0; i < nLocalPoints; ++i) {
//      std::cout << geomShapeFun0[i] << " "
//                << geomShapeFun1[i] << " "
//                << geomShapeFun2[i] << std::endl;
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
                           globalPoints.begin()+nEls*nLocalPoints,
                           globalPoints.begin()+2*nEls*nLocalPoints)),
      thrust::make_zip_iterator(
        thrust::make_tuple(globalPoints.begin()+nEls*nLocalPoints,
                           globalPoints.begin()+2*nEls*nLocalPoints,
                           globalPoints.end())),
      local2globalFunctor(nEls,
                          d_vtx0x.data(), d_vtx0y.data(), d_vtx0z.data(),
                          d_vtx1x.data(), d_vtx1y.data(), d_vtx1z.data(),
                          d_vtx2x.data(), d_vtx2y.data(), d_vtx2z.data(),
                          geomShapeFun0.data(), geomShapeFun1.data(), geomShapeFun2.data()));

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeMapping;
    cudaEventElapsedTime(&elapsedTimeMapping , start, stop);
    std::cout << "Time for mapping local to global coordinates is "
      << elapsedTimeMapping << " ms" << std::endl;

//    std::cout << "globalPoints = " << std::endl;
//    for (int i = 0; i < nLocalPoints; ++i) {
//      for (int j = 0; j < nEls; ++j) {
//        for (int k = 0; k < dim; ++k) {
//          std::cout << globalPoints[k * nLocalPoints * nEls + i * nEls + j] << " " << std::flush;
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
  }

} // namespace Bempp
