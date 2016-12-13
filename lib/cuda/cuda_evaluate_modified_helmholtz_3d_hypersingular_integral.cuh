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

#ifndef fiber_cuda_evaluate_modified_helmholtz_3d_hypersingular_integral_cuh
#define fiber_cuda_evaluate_modified_helmholtz_3d_hypersingular_integral_cuh

#include "cuda.cuh"
#include "cuda_evaluate_modified_helmholtz_3d_single_layer_potential_integral.cuh"

#include "../common/scalar_traits.hpp"

#include <device_launch_parameters.h>

namespace Fiber {

template<typename ValueType>
__device__ void Helmholtz3dHypersingularKernel(
    typename ScalarTraits<ValueType>::RealType testPointCoo[3],
    typename ScalarTraits<ValueType>::RealType trialPointCoo[3],
    const ValueType waveNumberImag,
    ValueType &result0Real, ValueType &result0Imag,
    ValueType &result1Real, ValueType &result1Imag) {

  Helmholtz3dSingleLayerPotentialKernel(
      testPointCoo, trialPointCoo, waveNumberImag,
      result0Real, result0Imag);

  ValueType factor = -(waveNumberImag * waveNumberImag);
  result1Real = result0Real * factor;
  result1Imag = result0Imag * factor;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateHelmholtz3dHypersingularIntegralFunctorCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int testIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testDofCount,
    const BasisFunctionType* __restrict__ testBasisValues,
    const BasisFunctionType* __restrict__ testBasisDerivatives,
    const unsigned int trialDofCount,
    const BasisFunctionType* __restrict__ trialBasisValues,
    const BasisFunctionType* __restrict__ trialBasisDerivatives,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* testElemNormals,
    const typename ScalarTraits<BasisFunctionType>::RealType* testIntegrationElements,
    const typename ScalarTraits<BasisFunctionType>::RealType* testJacobianInversesTransposed,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialElemNormals,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialIntegrationElements,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialJacobianInversesTransposed,
    const KernelType waveNumberImag, ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount) {

    const int coordCount = 3;

    // Determine trial and test element indices
    const int elemPairIndex = elemPairIndexBegin + i;
    const int trialElemPosition = trialIndices[elemPairIndex / testIndexCount];
    const int testElemPosition = testIndices[elemPairIndex % testIndexCount];

    // Gather integration elements
    const CoordinateType trialIntegrationElement =
        trialIntegrationElements[trialElemPosition];
    const CoordinateType testIntegrationElement =
        testIntegrationElements[testElemPosition];

    // Get normals
    CoordinateType trialElemNormal[coordCount], testElemNormal[coordCount];
#pragma unroll
    for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      trialElemNormal[coordIndex] =
          trialElemNormals[coordIndex * trialElemCount + trialElemPosition];
      testElemNormal[coordIndex] =
          testElemNormals[coordIndex * testElemCount + testElemPosition];
    }

    // Get transposed Jacobian inverses
    CoordinateType trialJacobianInverseTransposed[6], testJacobianInverseTransposed[6];
#pragma unroll
    for (int j = 0; j < 6; ++j) {
      trialJacobianInverseTransposed[j] =
          trialJacobianInversesTransposed[j * trialElemCount + trialElemPosition];
      testJacobianInverseTransposed[j] =
          testJacobianInversesTransposed[j * testElemCount + testElemPosition];
    }

    // Evaluate kernel
    KernelType kernelValues0Real[6 * 6], kernelValues0Imag[6 * 6],
        kernelValues1Real[6 * 6], kernelValues1Imag[6 * 6];
    CoordinateType trialPointCoo[coordCount], testPointCoo[coordCount];
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
        KernelType kernelValue0Real, kernelValue0Imag,
            kernelValue1Real, kernelValue1Imag;
#pragma unroll
        for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
          trialPointCoo[coordIndex] =
              trialGlobalPoints[coordIndex * trialPointCount * trialElemCount
                            + trialPoint * trialElemCount
                            + trialElemPosition];
          testPointCoo[coordIndex] =
              testGlobalPoints[coordIndex * testPointCount * testElemCount
                           + testPoint * testElemCount
                           + testElemPosition];
        }
        Helmholtz3dHypersingularKernel(
            testPointCoo, trialPointCoo, waveNumberImag,
            kernelValue0Real, kernelValue0Imag,
            kernelValue1Real, kernelValue1Imag);
        const size_t index = trialPoint * testPointCount + testPoint;
        kernelValues0Real[index] = kernelValue0Real;
        kernelValues0Imag[index] = kernelValue0Imag;
        kernelValues1Real[index] = kernelValue1Real;
        kernelValues1Imag[index] = kernelValue1Imag;
      }
    }

    // Perform numerical integration
    BasisFunctionType testSurfaceCurls[3], trialSurfaceCurls[3];
    BasisFunctionType vec[3];
    const size_t offsetResultImag = elemPairCount * trialDofCount * testDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        ResultType sumReal = 0., sumImag = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * constTrialQuadWeights[trialPoint];

          // Get surface curls
          BasisFunctionType trialBasisDerivativeDir0 =
              trialBasisDerivatives[0 * trialDofCount * trialPointCount
                                    + trialDof * trialPointCount
                                    + trialPoint];
          BasisFunctionType trialBasisDerivativeDir1 =
              trialBasisDerivatives[1 * trialDofCount * trialPointCount
                                    + trialDof * trialPointCount
                                    + trialPoint];
          vec[0] =
              trialBasisDerivativeDir0 * trialJacobianInverseTransposed[0] +
              trialBasisDerivativeDir1 * trialJacobianInverseTransposed[1];
          vec[1] =
              trialBasisDerivativeDir0 * trialJacobianInverseTransposed[2] +
              trialBasisDerivativeDir1 * trialJacobianInverseTransposed[3];
          vec[2] =
              trialBasisDerivativeDir0 * trialJacobianInverseTransposed[4] +
              trialBasisDerivativeDir1 * trialJacobianInverseTransposed[5];

          trialSurfaceCurls[0] = trialElemNormal[1] * vec[2] - trialElemNormal[2] * vec[1];
          trialSurfaceCurls[1] = trialElemNormal[2] * vec[0] - trialElemNormal[0] * vec[2];
          trialSurfaceCurls[2] = trialElemNormal[0] * vec[1] - trialElemNormal[1] * vec[0];

          ResultType partialSumReal = 0., partialSumImag = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[testPoint];

            // Get surface curls
            BasisFunctionType testBasisDerivativeDir0 =
                testBasisDerivatives[0 * testDofCount * testPointCount
                                     + testDof * testPointCount
                                     + testPoint];
            BasisFunctionType testBasisDerivativeDir1 =
                testBasisDerivatives[1 * testDofCount * testPointCount
                                     + testDof * testPointCount
                                     + testPoint];
            vec[0] =
                testBasisDerivativeDir0 * testJacobianInverseTransposed[0] +
                testBasisDerivativeDir1 * testJacobianInverseTransposed[1];
            vec[1] =
                testBasisDerivativeDir0 * testJacobianInverseTransposed[2] +
                testBasisDerivativeDir1 * testJacobianInverseTransposed[3];
            vec[2] =
                testBasisDerivativeDir0 * testJacobianInverseTransposed[4] +
                testBasisDerivativeDir1 * testJacobianInverseTransposed[5];

            testSurfaceCurls[0] = testElemNormal[1] * vec[2] - testElemNormal[2] * vec[1];
            testSurfaceCurls[1] = testElemNormal[2] * vec[0] - testElemNormal[0] * vec[2];
            testSurfaceCurls[2] = testElemNormal[0] * vec[1] - testElemNormal[1] * vec[0];

            const size_t indexKernelValues = trialPoint * testPointCount + testPoint;
            BasisFunctionType dotProduct0 = 0.;
            CoordinateType dotProduct1 = 0.;
#pragma unroll
            for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
              dotProduct0 +=
                  testSurfaceCurls[coordIndex] * trialSurfaceCurls[coordIndex];
              dotProduct1 +=
                  testElemNormal[coordIndex] * trialElemNormal[coordIndex];
            }
            ResultType term0Real =
                kernelValues0Real[indexKernelValues] * dotProduct0;
            ResultType term0Imag =
                kernelValues0Imag[indexKernelValues] * dotProduct0;

            ResultType term1Real =
                kernelValues1Real[indexKernelValues] * dotProduct1
                * testBasisValues[testDof * testPointCount + testPoint]
                * trialBasisValues[trialDof * trialPointCount + trialPoint];
            ResultType term1Imag =
                kernelValues1Imag[indexKernelValues] * dotProduct1
                * testBasisValues[testDof * testPointCount + testPoint]
                * trialBasisValues[trialDof * trialPointCount + trialPoint];

            ResultType integrandReal = term0Real + term1Real;
            ResultType integrandImag = term0Imag + term1Imag;

            partialSumReal += integrandReal * testWeight;
            partialSumImag += integrandImag * testWeight;
          }
          sumReal += partialSumReal * trialWeight;
          sumImag += partialSumImag * trialWeight;
        }
        const size_t indexResult = trialDof * testDofCount * elemPairCount
                                   + testDof * elemPairCount
                                   + i;
        result[indexResult] = sumReal;
        result[indexResult + offsetResultImag] = sumImag;
      }
    }
  }
}

} // namespace Fiber

#endif
