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

#ifndef fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_cuh
#define fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_cuh

//#define CUDA_QUAD_WEIGHTS {1.11690794839005499983e-01, 5.49758718276609978370e-02, 1.11690794839005499983e-01, 1.11690794839005499983e-01, 5.49758718276609978370e-02, 5.49758718276609978370e-02}
//#define CUDA_BASIS_VALUES {1.08103014528751373291e-01, 8.16847562789916992188e-01, 4.45948481559753417969e-01, 4.45948481559753417969e-01, 9.15762111544609069824e-02, 9.15762111544609069824e-02, 4.45948481559753417969e-01, 9.15762111544609069824e-02, 1.08103014528751373291e-01, 4.45948481559753417969e-01, 8.16847562789916992188e-01, 9.15762111544609069824e-02, 4.45948481559753417969e-01, 9.15762111544609069824e-02, 4.45948481559753417969e-01, 1.08103014528751373291e-01, 9.15762111544609069824e-02, 8.16847562789916992188e-01}

#include "../common/scalar_traits.hpp"

#include <device_launch_parameters.h>

__constant__ double constTestQuadWeights[3][6];
__constant__ double constTrialQuadWeights[3][6];

__constant__ double constTestGeomShapeFun0[6];
__constant__ double constTestGeomShapeFun1[6];
__constant__ double constTestGeomShapeFun2[6];

__constant__ double constTrialGeomShapeFun0[6];
__constant__ double constTrialGeomShapeFun1[6];
__constant__ double constTrialGeomShapeFun2[6];

namespace Fiber {

template<typename ValueType>
inline __device__ void Laplace3dSingleLayerPotentialKernel(
    typename ScalarTraits<ValueType>::RealType testPointCoo[3],
    typename ScalarTraits<ValueType>::RealType trialPointCoo[3],
    ValueType* result) {

  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  static constexpr int COORD_COUNT = 3;

  ValueType sum = 0;
  for (int coordIndex = 0; coordIndex < COORD_COUNT; ++coordIndex) {
    ValueType diff = testPointCoo[coordIndex] - trialPointCoo[coordIndex];
    sum += diff * diff;
  }
  result[0] = 0.25 * static_cast<CoordinateType>(M_1_PI) * rsqrt(sum);

}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          int testPointCount_, int trialPointCount_,
          int testDofCount_  , int trialDofCount_>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int testIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const BasisFunctionType* __restrict__ testBasisValues,
    const BasisFunctionType* __restrict__ trialBasisValues,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* testIntegrationElements,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialIntegrationElements,
    ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  static constexpr int       COORD_COUNT =                3;
  static constexpr int  TEST_POINT_COUNT =  testPointCount_;
  static constexpr int TRIAL_POINT_COUNT = trialPointCount_;
  static constexpr int    TEST_DOF_COUNT =    testDofCount_;
  static constexpr int   TRIAL_DOF_COUNT =   trialDofCount_;

//  constexpr ResultType        QUAD_WEIGHTS[  6] = CUDA_QUAD_WEIGHTS;
//  constexpr BasisFunctionType BASIS_VALUES[3*6] = CUDA_BASIS_VALUES;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount) {

    // Determine trial and test element indices
    const int     elemPairIndex = elemPairIndexBegin + i;
    const int trialElemPosition = trialIndices[elemPairIndex / testIndexCount];
    const int  testElemPosition =  testIndices[elemPairIndex % testIndexCount];

    // Gather integration elements
    const CoordinateType trialIntegrationElement =
        trialIntegrationElements[trialElemPosition];
    const CoordinateType  testIntegrationElement =
        testIntegrationElements [ testElemPosition];

    // Evaluate kernel
    KernelType kernelValues[TEST_POINT_COUNT * TRIAL_POINT_COUNT];
    CoordinateType trialPointCoo[COORD_COUNT], testPointCoo[COORD_COUNT];
    for (int trialPoint = 0; trialPoint < TRIAL_POINT_COUNT; ++trialPoint) {
      for (int testPoint = 0; testPoint < TEST_POINT_COUNT; ++testPoint) {
        KernelType kernelValue;
        for (int coordIndex = 0; coordIndex < COORD_COUNT; ++coordIndex) {
          trialPointCoo[coordIndex] =
              trialGlobalPoints[coordIndex * TRIAL_POINT_COUNT * trialElemCount
                              + trialPoint * trialElemCount
                              + trialElemPosition];
          testPointCoo[coordIndex] =
              testGlobalPoints[coordIndex * TEST_POINT_COUNT * testElemCount
                             + testPoint * testElemCount
                             + testElemPosition];
        }
        Laplace3dSingleLayerPotentialKernel(testPointCoo, trialPointCoo,
            &kernelValue);
        kernelValues[trialPoint * TEST_POINT_COUNT + testPoint] = kernelValue;
      }
    }

    // Perform numerical integration
    for (int trialDof = 0; trialDof < TRIAL_DOF_COUNT; ++trialDof) {
      for (int testDof = 0; testDof < TEST_DOF_COUNT; ++testDof) {
        ResultType sum = 0.;
        for (int trialPoint = 0; trialPoint < TRIAL_POINT_COUNT; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * constTrialQuadWeights[0][trialPoint];
          ResultType partialSum = 0.;
          for (int testPoint = 0; testPoint < TEST_POINT_COUNT; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[0][testPoint];
            partialSum += kernelValues[trialPoint *  TEST_POINT_COUNT +  testPoint]
                    *  testBasisValues[   testDof *  TEST_POINT_COUNT +  testPoint]
                    * trialBasisValues[  trialDof * TRIAL_POINT_COUNT + trialPoint]
                    * testWeight;
          }
          sum += partialSum * trialWeight;
        }

        // Write result to global memory
        const size_t index = trialDof * TEST_DOF_COUNT * elemPairCount
                           +  testDof * elemPairCount
                           + i;
        result[index] = sum;
      }
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          int testPointCount_, int trialPointCount_>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialDofFunctor(
    const int*  __restrict__ testElemIndices,
    const char* __restrict__ testLocalDofIndices ,
    const size_t             numberOfTestDofs,
    const int*  __restrict__ trialElemIndices,
    const char* __restrict__ trialLocalDofIndices,
    const size_t             numberOfTrialDofs,
    const BasisFunctionType* __restrict__ testBasisValues,
    const BasisFunctionType* __restrict__ trialBasisValues,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* testIntegrationElements,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialIntegrationElements,
    ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  static constexpr int       COORD_COUNT =                3;
  static constexpr int  TEST_POINT_COUNT =  testPointCount_;
  static constexpr int TRIAL_POINT_COUNT = trialPointCount_;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < numberOfTestDofs * numberOfTrialDofs) {

    // Determine pair of trial and test dofs to compute
    const int trialDof = i / numberOfTestDofs;
    const int  testDof = i % numberOfTestDofs;

    ResultType coeff = ResultType(0);

    for (int trialElem = 0; trialElem < 10; ++trialElem) {

      const int trialElemIdx     = trialElemIndices    [trialDof * 10 + trialElem];
      if (trialElemIdx == -1) break;
      const int trialLocalDofIdx = trialLocalDofIndices[trialDof * 10 + trialElem];
      const CoordinateType trialIntegrationElement =
          trialIntegrationElements[trialElemIdx];

      for (int testElem = 0; testElem < 10; ++testElem) {

        const int testElemIdx = testElemIndices[testDof * 10 + testElem];
        if (testElemIdx == -1) break;
        const int testLocalDofIdx = testLocalDofIndices[testDof * 10 + testElem];
        const CoordinateType testIntegrationElement =
            testIntegrationElements[testElemIdx];

        // Evaluate kernel
        KernelType kernelValues[TEST_POINT_COUNT * TRIAL_POINT_COUNT];
        CoordinateType trialPointCoo[COORD_COUNT], testPointCoo[COORD_COUNT];
        for (int trialPoint = 0; trialPoint < TRIAL_POINT_COUNT; ++trialPoint) {
          for (int testPoint = 0; testPoint < TEST_POINT_COUNT; ++testPoint) {
            KernelType kernelValue;
            for (int coordIndex = 0; coordIndex < COORD_COUNT; ++coordIndex) {
              trialPointCoo[coordIndex] =
                  trialGlobalPoints[coordIndex * TRIAL_POINT_COUNT * trialElemCount
                                  + trialPoint * trialElemCount
                                  + trialElemIdx];
              testPointCoo[coordIndex] =
                  testGlobalPoints[coordIndex * TEST_POINT_COUNT * testElemCount
                                 + testPoint * testElemCount
                                 + testElemIdx];
            }
            Laplace3dSingleLayerPotentialKernel(testPointCoo, trialPointCoo,
                &kernelValue);
            kernelValues[trialPoint * TEST_POINT_COUNT + testPoint] = kernelValue;
          }
        }

        ResultType sum = ResultType(0);
        for (int trialPoint = 0; trialPoint < TRIAL_POINT_COUNT; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * constTrialQuadWeights[0][trialPoint];
          ResultType partialSum = 0.;
          for (int testPoint = 0; testPoint < TEST_POINT_COUNT; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[0][testPoint];
            partialSum +=
                kernelValues[trialPoint * TEST_POINT_COUNT + testPoint]
                *  testBasisValues[ testLocalDofIdx *  TEST_POINT_COUNT +  testPoint]
                * trialBasisValues[trialLocalDofIdx * TRIAL_POINT_COUNT + trialPoint]
                * testWeight;
          } // for testPoint
          sum += partialSum * trialWeight;
        } // for trialPoint

        coeff += sum;

      } // for testElem
    } // for trialElem

    result[i] = coeff;
  }
}

} // namespace Fiber

#endif
