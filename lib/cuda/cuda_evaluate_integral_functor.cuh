// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_cuda_evaluate_integral_functor_cuh
#define fiber_cuda_evaluate_integral_functor_cuh

#include "cuda.cuh"

#include "../common/common.hpp"
#include "../common/scalar_traits.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct CudaEvaluateIntegralFunctorCached {

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

protected:
  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData<CoordinateType> testQuadData, trialQuadData;
  BasisFunData<BasisFunctionType> testBasisData, trialBasisData;
  thrust::device_ptr<ResultType> result;

  ElemData<CoordinateType> testElemData, trialElemData;

public:
  CudaEvaluateIntegralFunctorCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData<CoordinateType> _testQuadData,
      const QuadData<CoordinateType> _trialQuadData,
      const BasisFunData<BasisFunctionType> _testBasisData,
      const BasisFunData<BasisFunctionType> _trialBasisData,
      const ElemData<CoordinateType> _testElemData,
      const ElemData<CoordinateType> _trialElemData,
      thrust::device_ptr<ResultType> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testElemData(_testElemData), trialElemData(_trialElemData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testBasisData.dofCount * trialBasisData.dofCount <= 36);
  }

  __device__
  virtual void operator() (const unsigned int i) = 0;
};

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct CudaEvaluateIntegralFunctorNonCached {

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

protected:
  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData<CoordinateType> testQuadData, trialQuadData;
  BasisFunData<BasisFunctionType> testBasisData, trialBasisData;
  thrust::device_ptr<ResultType> result;

  RawGeometryData<double> testRawGeometryData, trialRawGeometryData;
  GeomShapeFunData<CoordinateType> testGeomShapeFunData, trialGeomShapeFunData;

public:
  CudaEvaluateIntegralFunctorNonCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData<CoordinateType> _testQuadData,
      const QuadData<CoordinateType> _trialQuadData,
      const BasisFunData<BasisFunctionType> _testBasisData,
      const BasisFunData<BasisFunctionType> _trialBasisData,
      const RawGeometryData<double> _testRawGeometryData,
      const RawGeometryData<double> _trialRawGeometryData,
      const GeomShapeFunData<CoordinateType> _testGeomShapeFunData,
      const GeomShapeFunData<CoordinateType> _trialGeomShapeFunData,
      thrust::device_ptr<ResultType> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testRawGeometryData(_testRawGeometryData),
        trialRawGeometryData(_trialRawGeometryData),
        testGeomShapeFunData(_testGeomShapeFunData),
        trialGeomShapeFunData(_trialGeomShapeFunData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testBasisData.dofCount * trialBasisData.dofCount <= 36);
  }

  __device__
  virtual void operator() (const unsigned int i) = 0;
};

} // namespace Fiber

#endif
