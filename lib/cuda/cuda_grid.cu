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
#include <thrust/iterator/counting_iterator.h>

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

    unsigned int elemCount;

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
      const unsigned int _elemCount,
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
      : elemCount(_elemCount),
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
      unsigned int localPointIdx = i / elemCount;
      unsigned int elementIdx = i % elemCount;

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

  struct evaluateKernelFunctor {

    unsigned int testPointCount;
    unsigned int trialPointCount;

    thrust::device_ptr<double> testPoints;
    thrust::device_ptr<double> trialPoints;

    evaluateKernelFunctor(
      const unsigned int _testPointCount, const unsigned int _trialPointCount,
      thrust::device_ptr<double> _testPoints,
      thrust::device_ptr<double> _trialPoints)
      : testPointCount(_testPointCount), trialPointCount(_trialPointCount),
        testPoints(_testPoints), trialPoints(_trialPoints) {}

    __host__ __device__
    double operator() (const unsigned int i) {

      // Which memory mapping to go for?
//      unsigned int testPointIdx = i % nTestPoints;
//      unsigned int trialPointIdx = i / nTestPoints;
      unsigned int testPointIdx = i / trialPointCount;
      unsigned int trialPointIdx = i % trialPointCount;

      double xTest = testPoints[testPointIdx];
      double yTest = testPoints[testPointIdx+testPointCount];
      double zTest = testPoints[testPointIdx+2*testPointCount];

      double xTrial = trialPoints[trialPointIdx];
      double yTrial = trialPoints[trialPointIdx+trialPointCount];
      double zTrial = trialPoints[trialPointIdx+2*trialPointCount];

      double xDiff = xTest - xTrial;
      double yDiff = yTest - yTrial;
      double zDiff = zTest - zTrial;

      double distance = std::sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);

      return 1.0 / (4.0 * M_PI * distance);
    }
  };

  struct evaluateIntegralFunctor {

    unsigned int testElemCount;
    unsigned int trialElemCount;
    unsigned int localTestPointCount;
    unsigned int localTrialPointCount;
    unsigned int localTestDofCount;
    unsigned int localTrialDofCount;

    thrust::device_ptr<double> kernelValues;
    thrust::device_ptr<double> integrationElements;
    thrust::device_ptr<double> testBasisData;
    thrust::device_ptr<double> trialBasisData;
    thrust::device_ptr<double> testQuadWeights;
    thrust::device_ptr<double> trialQuadWeights;
    thrust::device_ptr<double> result;

    evaluateIntegralFunctor(
      const unsigned int _testElemCount, const unsigned int _trialElemCount,
      const unsigned int _localTestPointCount,
      const unsigned int _localTrialPointCount,
      const unsigned int _localTestDofCount,
      const unsigned int _localTrialDofCount,
      thrust::device_ptr<double> _kernelValues,
      thrust::device_ptr<double> _integrationElements,
      thrust::device_ptr<double> _testBasisData,
      thrust::device_ptr<double> _trialBasisData,
      thrust::device_ptr<double> _testQuadWeights,
      thrust::device_ptr<double> _trialQuadWeights,
      thrust::device_ptr<double> _result)
      : testElemCount(_testElemCount), trialElemCount(_trialElemCount),
        localTestPointCount(_localTestPointCount),
        localTrialPointCount(_localTrialPointCount),
        localTestDofCount(_localTestDofCount),
        localTrialDofCount(_localTrialDofCount),
        kernelValues(_kernelValues),
        integrationElements(_integrationElements),
        testBasisData(_testBasisData),
        trialBasisData(_trialBasisData),
        testQuadWeights(_testQuadWeights),
        trialQuadWeights(_trialQuadWeights),
        result(_result) {}

    __host__ __device__
    void operator() (const unsigned int i) {

      // Identify element pair to work on
      const unsigned int testElementIdx = i / trialElemCount;
      const unsigned int trialElementIdx = i % trialElemCount;

      // Fetch element pair related data
      double testIntegrationElement = integrationElements[testElementIdx];
      double trialIntegrationElement = integrationElements[trialElementIdx];
      double* elementPairKernelValues =
          new double[localTestPointCount * localTrialPointCount];
      for (int tePt = 0; tePt < localTestPointCount; ++tePt) {
        for (int trPt = 0; trPt < localTrialPointCount; ++trPt) {
          unsigned int localPointPairIdx = tePt * localTestPointCount + trPt;
          unsigned int globalPointPairIdx =
            testElementIdx * trialElemCount * localTestPointCount * localTrialPointCount
            + trialElementIdx * localTrialPointCount
            + tePt * localTrialPointCount * trialElemCount
            + trPt;
          elementPairKernelValues[localPointPairIdx] =
              kernelValues[globalPointPairIdx];
        }
      }

      // Create element pair result array
      double* elementPairResult =
          new double[localTestDofCount * localTrialDofCount];

      // Loop over dofs
      for (int teDof = 0; teDof < localTestDofCount; ++teDof) {
        for (int trDof = 0; trDof < localTrialDofCount; ++trDof) {
          double sum = 0.0;
          // Loop over quad points
          for (int tePt = 0; tePt < localTestPointCount; ++tePt) {
            double testWeight = testIntegrationElement * testQuadWeights[tePt];
            double partialSum = 0.0;
            for (int trPt = 0; trPt < localTrialPointCount; ++trPt) {
              double trialWeight = trialIntegrationElement * testQuadWeights[trPt];
              partialSum +=
                  trialWeight
                  * testBasisData[teDof * localTestDofCount + tePt]
                  * trialBasisData[trDof * localTrialDofCount + trPt]
                  * elementPairKernelValues[tePt * localTestPointCount + trPt];
            }
            sum += partialSum * testWeight;
          }
          elementPairResult[teDof * localTestDofCount + trDof] = sum;
        }
      }

      // Copy data to global result array
      unsigned int offset = i * localTestDofCount * localTrialDofCount;
      for (int teDof = 0; teDof < localTestDofCount; ++teDof) {
        for (int trDof = 0; trDof < localTrialDofCount; ++trDof) {
          unsigned int dofPairIdx = teDof * localTestDofCount + trDof;
          result[offset + dofPairIdx] = elementPairResult[dofPairIdx];
        }
      }
      delete[] elementPairKernelValues;
      delete[] elementPairResult;
    }
  };

  CudaGrid::CudaGrid() {

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

    m_normals.clear();
    m_integrationElements.clear();

    m_setupDone = false;
  }

  CudaGrid::~CudaGrid() {

  }

  void CudaGrid::pushGeometry(const Matrix<double> &vertices,
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

    m_setupDone = false;
  }

  void CudaGrid::setupGeometry() {

    if (m_setupDone == false) {

      // Gather element corner coordinates
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

  //    std::cout << "m_vtx0x = " << std::endl;
  //    for (int i = 0; i < m_nEls; ++i) {
  //      std::cout << m_vtx0x[i] << std::endl;
  //    }
  //    std::cout << std::endl;

      // TODO Remove m_vertices and m_elementCorners on the device at this point?

      calculateNormalsAndIntegrationElements();

      m_setupDone = true;
    }

    // TEST Convert local to global coordinates for all elements
    thrust::host_vector<double> h_localPoints(2);
    h_localPoints[0] = 0.1;
    h_localPoints[1] = 0.1;

    thrust::device_vector<double> d_globalPoints;
    local2global(h_localPoints, d_globalPoints);

    // TEST Evaluate the kernel
    thrust::device_vector<double> d_kernelValues;
    evaluateKernel(d_globalPoints, d_globalPoints, d_kernelValues);

    // TEST Evaluate the integral
    thrust::device_vector<double> d_testBasisData(3);
    d_testBasisData[0] = 0.9;
    d_testBasisData[1] = 0.2;
    d_testBasisData[2] = 0.2;
    thrust::device_vector<double> d_trialBasisData = d_testBasisData;

    thrust::device_vector<double> d_testQuadWeights(1);
    d_testQuadWeights[0] = 2.0;
    thrust::device_vector<double> d_trialQuadWeights = d_testQuadWeights;

    thrust::device_vector<double> d_result;
    evaluateIntegral(m_ElemCount, m_ElemCount, d_kernelValues,
        d_testBasisData, d_trialBasisData, d_testQuadWeights, d_trialQuadWeights, d_result);
  }

  void CudaGrid::calculateNormalsAndIntegrationElements() {

    // Allocate device memory
    m_normals.resize(m_dim * m_ElemCount);
    m_integrationElements.resize(m_ElemCount);

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
        thrust::make_tuple(m_normals.begin(), m_normals.begin()+m_ElemCount,
                           m_normals.begin()+2*m_ElemCount, m_integrationElements.begin())),
      calculateElementNormalAndIntegrationElementFunctor());

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeNormals;
    cudaEventElapsedTime(&elapsedTimeNormals , start, stop);
    std::cout << "Time for calculating normals and integration elements is "
      << elapsedTimeNormals << " ms" << std::endl;

//    std::cout << "m_normals = " << std::endl;
//    for (int i = 0; i < m_nEls; ++i) {
//      for (int j = 0; j < m_dim; ++j) {
//        std::cout << m_normals[j * m_nEls + i] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    std::cout << "m_integrationElements = " << std::endl;
//    for (int i = 0; i < m_nEls; ++i) {
//      std::cout << m_integrationElements[i] << std::endl;
//    }
//    std::cout << std::endl;
  }

  void CudaGrid::local2global(const thrust::host_vector<double> &h_localPoints,
                              thrust::device_vector<double> &d_globalPoints) {

    // localPoints = [r0 r1 ... rN | s0 s1 ... sN]

    const unsigned int localPointDim = 2;
    const unsigned int localPointCount = h_localPoints.size() / localPointDim;

    if (h_localPoints.size() % localPointDim != 0)
      throw std::runtime_error("CudaGrid::local2global(): "
                               "only valid for two-dimensional local points");

    // Allocate device memory
    d_globalPoints.resize(m_dim * m_ElemCount * localPointCount);
    // [xPt0el0 xPt0el1 ... xPt0elM | xPt1el0 ... xPt1elM | ... | ... xPtNelM |
    //  yPt0el0 yPt0el1 ... yPt0elM | yPt1el0 ... yPt1elM | ... | ... yPtNelM |
    //  zPt0el0 zPt0el1 ... zPt0elM | zPt1el0 ... zPt1elM | ... | ... zPtNelM ]

    thrust::host_vector<double> h_geomShapeFun0(localPointCount);
    thrust::host_vector<double> h_geomShapeFun1(localPointCount);
    thrust::host_vector<double> h_geomShapeFun2(localPointCount);
    // [pt0 pt1 ... ptN]

    // Evaluate geometrical shape function values on the host
    thrust::transform(thrust::host,
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localPoints.begin(),
                           h_localPoints.begin()+localPointCount)),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localPoints.begin()+localPointCount,
                           h_localPoints.end())),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_geomShapeFun0.begin(),
                           h_geomShapeFun1.begin(),
                           h_geomShapeFun2.begin())),
      geomShapeFunFunctor());

    // Copy data to device
    thrust::device_vector<double> d_geomShapeFun0 = h_geomShapeFun0;
    thrust::device_vector<double> d_geomShapeFun1 = h_geomShapeFun1;
    thrust::device_vector<double> d_geomShapeFun2 = h_geomShapeFun2;

//    std::cout << "d_geomShapeFun = " << std::endl;
//    for (int i = 0; i < nLocalPoints; ++i) {
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
        thrust::make_tuple(d_globalPoints.begin(),
                           d_globalPoints.begin()+m_ElemCount*localPointCount,
                           d_globalPoints.begin()+2*m_ElemCount*localPointCount)),
      thrust::make_zip_iterator(
        thrust::make_tuple(d_globalPoints.begin()+m_ElemCount*localPointCount,
                           d_globalPoints.begin()+2*m_ElemCount*localPointCount,
                           d_globalPoints.end())),
      local2globalFunctor(m_ElemCount,
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

//    std::cout << "d_globalPoints = " << std::endl;
//    for (int i = 0; i < nLocalPoints; ++i) {
//      for (int j = 0; j < m_nEls; ++j) {
//        for (int k = 0; k < m_dim; ++k) {
//          std::cout << d_globalPoints[k * localPointsCount * m_nEls + i * m_nEls + j] << " " << std::flush;
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
  }

  void CudaGrid::evaluateKernel(thrust::device_vector<double> &d_testPoints,
                                thrust::device_vector<double> &d_trialPoints,
                                thrust::device_vector<double> &d_kernelValues) {

    const unsigned int pointDim = 3;
    const unsigned int testPointCount = d_testPoints.size() / pointDim;
    const unsigned int trialPointCount = d_trialPoints.size() / pointDim;

    if (d_testPoints.size() % pointDim != 0 || d_trialPoints.size() % pointDim != 0)
      throw std::runtime_error("CudaGrid::evaluateKernel(): "
                               "only valid for three-dimensional points");

    // Allocate device memory
    d_kernelValues.resize(testPointCount * trialPointCount);

    // Measure time of the GPU execution (CUDA event based)
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    thrust::tabulate(d_kernelValues.begin(), d_kernelValues.end(),
                     evaluateKernelFunctor(testPointCount, trialPointCount,
                         d_testPoints.data(), d_trialPoints.data()));

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeKernel;
    cudaEventElapsedTime(&elapsedTimeKernel , start, stop);
    std::cout << "Time for kernel evaluation is "
      << elapsedTimeKernel << " ms" << std::endl;

//    std::cout << "d_kernelValues = " << std::endl;
//    for (int i = 0; i < nTestPoints; ++i) {
//      for (int j = 0; j < nTrialPoints; ++j) {
//        std::cout << d_kernelValues[i * nTestPoints + j] << " " << std::flush;
//      }
//      std::cout << std::endl;
//    }
  }

  void CudaGrid::evaluateIntegral(
      const unsigned int testElemCount, const unsigned int trialElemCount,
      thrust::device_vector<double> &d_kernelValues,
      thrust::device_vector<double> &d_testBasisData,
      thrust::device_vector<double> &d_trialBasisData,
      thrust::device_vector<double> &d_testQuadWeights,
      thrust::device_vector<double> &d_trialQuadWeights,
      thrust::device_vector<double> &d_result) {

    const unsigned int localTestPointCount = d_testQuadWeights.size();
    const unsigned int localTrialPointCount = d_trialQuadWeights.size();

    const unsigned int localTestDofCount =
        d_testBasisData.size() / localTestPointCount;
    const unsigned int localTrialDofCount =
        d_trialBasisData.size() / localTrialPointCount;

    if (d_testBasisData.size() % localTestPointCount != 0
        || d_trialBasisData.size() % localTrialPointCount != 0)
      throw std::runtime_error("CudaGrid::evaluateIntegral(): "
                               "basis data size does not fit number of quad points");

    // Allocate device memory
    d_result.resize(testElemCount * localTestDofCount
        * trialElemCount * localTrialDofCount);

    // Measure time of the GPU execution (CUDA event based)
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+testElemCount*trialElemCount,
                     evaluateIntegralFunctor(testElemCount, trialElemCount,
                                             localTestPointCount, localTrialPointCount,
                                             localTestDofCount, localTrialDofCount,
                                             m_integrationElements.data(),
                                             d_kernelValues.data(),
                                             d_testBasisData.data(),
                                             d_trialBasisData.data(),
                                             d_testQuadWeights.data(),
                                             d_trialQuadWeights.data(),
                                             d_result.data()));

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float elapsedTimeIntegral;
    cudaEventElapsedTime(&elapsedTimeIntegral , start, stop);
    std::cout << "Time for integral evaluation is "
      << elapsedTimeIntegral << " ms" << std::endl;

//    std::cout << "d_result = " << std::endl;
//    for (int i = 0; i < nTestElements; ++i) {
//      for (int j = 0; j < nTrialElements; ++j) {
//        for (int k = 0; k < nLocalTestDofs; ++k) {
//          for (int l = 0; l < nLocalTrialDofs; ++l) {
//            std::cout << d_result[
//              i * nTrialElements * nLocalTestDofs * nLocalTrialDofs
//              + j * nLocalTestDofs * nLocalTrialDofs
//              + k * nLocalTrialDofs
//              +l] << " " << std::flush;
//          }
//        }
//      }
//      std::cout << std::endl;
//    }
  }

} // namespace Bempp
