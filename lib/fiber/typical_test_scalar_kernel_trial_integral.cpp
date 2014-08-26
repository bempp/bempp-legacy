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

#ifndef fiber_default_test_single_scalar_kernel_trial_integral_imp_hpp
#define fiber_default_test_single_scalar_kernel_trial_integral_imp_hpp

#include "typical_test_scalar_kernel_trial_integral.hpp"

#include "collection_of_4d_arrays.hpp"
#include "explicit_instantiation.hpp"
#include "geometrical_data.hpp"
#include "../common/acc.hpp"
#include "../common/complex_aux.hpp"

#include <cassert>
#include <iostream>
#include <tbb/scalable_allocator.h>

namespace Fiber {

// Internal functions
namespace {

template <typename BasisFunctionType, typename ResultType>
void setToProduct(const arma::Mat<BasisFunctionType> &source,
                  BasisFunctionType weight, arma::Mat<ResultType> &result) {
  result = weight * source.t();
}

template <typename BasisFunctionType, typename ResultType>
void setToProduct(const arma::Mat<BasisFunctionType> &source,
                  typename ScalarTraits<BasisFunctionType>::RealType weight,
                  arma::Mat<ResultType> &result) {
  result = weight * source.t();
}

template <typename CoordinateType>
void setToProduct(const arma::Mat<CoordinateType> &source,
                  CoordinateType weight,
                  arma::Mat<std::complex<CoordinateType>> &result) {
  result.fill(0.);
  result.set_real(weight * source.t());
}

template <typename BasisFunctionType, typename ResultType>
void outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
    const _3dArray<BasisFunctionType> &values,
    const GeometricalData<typename ScalarTraits<BasisFunctionType>::RealType> &
        geomData,
    const std::vector<typename ScalarTraits<BasisFunctionType>::RealType> &
        quadWeights,
    std::vector<ResultType, tbb::scalable_allocator<ResultType>> &tmp) {
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const size_t transDim = values.extent(0);
  const size_t dofCount = values.extent(1);
  const size_t pointCount = values.extent(2);
  const size_t sliceSize = transDim * dofCount;

  const BasisFunctionType *sourceSliceStart = values.begin();
  tmp.resize(dofCount * transDim * pointCount);
  ResultType *destSliceStart = &tmp[0];

  for (size_t point = 0; point < pointCount; ++point) {
    arma::Mat<BasisFunctionType> origMat(
        const_cast<BasisFunctionType *>(sourceSliceStart), transDim, dofCount,
        false);
    arma::Mat<ResultType> transposedMat(destSliceStart, dofCount, transDim,
                                        false);
    const CoordinateType weight =
        geomData.integrationElements(point) * quadWeights[point];
    // we take the complex conjugate here
    setToProduct(origMat, weight, transposedMat);
    sourceSliceStart += sliceSize;
    destSliceStart += sliceSize;
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
    const _3dArray<BasisFunctionType> &values,
    const std::vector<KernelType, tbb::scalable_allocator<KernelType>> &weights,
    std::vector<ResultType, tbb::scalable_allocator<ResultType>> &tmp) {
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const size_t transDim = values.extent(0);
  const size_t dofCount = values.extent(1);
  const size_t pointCount = values.extent(2);
  const size_t sliceSize = transDim * dofCount;

  const BasisFunctionType *sourceSliceStart = values.begin();
  tmp.resize(dofCount * transDim * pointCount);
  ResultType *destSliceStart = &tmp[0];

  for (size_t point = 0; point < pointCount; ++point) {
    arma::Mat<BasisFunctionType> origMat(
        const_cast<BasisFunctionType *>(sourceSliceStart), transDim, dofCount,
        false);
    arma::Mat<ResultType> transposedMat(destSliceStart, dofCount, transDim,
                                        false);
    const KernelType weight = weights[point];
    // we take the complex conjugate here
    setToProduct(origMat, weight, transposedMat);
    sourceSliceStart += sliceSize;
    destSliceStart += sliceSize;
  }
}

template <typename BasisFunctionType, typename ResultType>
void outOfPlaceConjugateTransposeDimAndDofDimensions(
    const _3dArray<BasisFunctionType> &values,
    std::vector<ResultType, tbb::scalable_allocator<ResultType>> &tmp) {
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const size_t transDim = values.extent(0);
  const size_t dofCount = values.extent(1);
  const size_t pointCount = values.extent(2);
  const size_t sliceSize = transDim * dofCount;

  const BasisFunctionType *sourceSliceStart = values.begin();
  tmp.resize(dofCount * transDim * pointCount);
  ResultType *destSliceStart = &tmp[0];

  for (size_t point = 0; point < pointCount; ++point) {
    arma::Mat<BasisFunctionType> origMat(
        const_cast<BasisFunctionType *>(sourceSliceStart), transDim, dofCount,
        false);
    arma::Mat<ResultType> transposedMat(destSliceStart, dofCount, transDim,
                                        false);
    // we take the complex conjugate here
    transposedMat = origMat.t();
    sourceSliceStart += sliceSize;
    destSliceStart += sliceSize;
  }
}

template <typename BasisFunctionType, typename ResultType>
inline void addABt(const arma::Mat<BasisFunctionType> &A,
                   const arma::Mat<BasisFunctionType> &B,
                   arma::Mat<ResultType> &C) {
  C = C + A * B.t();
}

template <typename ResultType>
inline void addABt(const arma::Mat<ResultType> &A,
                   const arma::Mat<ResultType> &B, arma::Mat<ResultType> &C) {
  C += A * B.t();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void evaluateWithNontensorQuadratureRuleStandardImpl(
    const GeometricalData<typename ScalarTraits<ResultType>::RealType> &
        testGeomData,
    const GeometricalData<typename ScalarTraits<ResultType>::RealType> &
        trialGeomData,
    const CollectionOf3dArrays<BasisFunctionType> &testValues,
    const CollectionOf3dArrays<BasisFunctionType> &trialValues,
    const CollectionOf3dArrays<KernelType> &kernelValues,
    const std::vector<typename ScalarTraits<ResultType>::RealType> &quadWeights,
    arma::Mat<ResultType> &result) {
  // We assume the integrand has the structure
  // (i -- transformations)
  // sum_i (\vec test_i \cdot \vec trial_i) kernel
  // or sum_i (\vec test_i \cdot kernel_i \vec trial_i) kernel

  // Evaluate constants and assert that array dimensions are correct
  const size_t transCount = testValues.size();
  assert(transCount >= 1);
  assert(trialValues.size() == transCount);
  assert(kernelValues.size() == 1 || kernelValues.size() == transCount);

  const size_t testDofCount = testValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(testValues[i].extent(1) == testDofCount);
  const size_t trialDofCount = trialValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(trialValues[i].extent(1) == trialDofCount);
  const size_t pointCount = quadWeights.size();

  for (size_t i = 0; i < kernelValues.size(); ++i)
    assert(kernelValues[i].extent(2) == pointCount);
  for (size_t i = 0; i < transCount; ++i)
    assert(testValues[i].extent(2) == pointCount);
  for (size_t i = 0; i < transCount; ++i)
    assert(trialValues[i].extent(2) == pointCount);
  assert(result.n_rows == testDofCount);
  assert(result.n_cols == trialDofCount);

  result.fill(0);

  std::vector<KernelType, tbb::scalable_allocator<KernelType>> products(
      pointCount);
  std::vector<BasisFunctionType, tbb::scalable_allocator<BasisFunctionType>>
  tmpTest, tmpTrial;

  for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
    const size_t transDim = testValues[transIndex].extent(0);
    assert(trialValues[transIndex].extent(0) == transDim);

    if (transIndex == 0 || kernelValues.size() > 1)
      for (size_t point = 0; point < pointCount; ++point)
        // we take the complex conj. here, later we'll remove it
        products[point] = conj(kernelValues[transIndex](0, 0, point)) *
                          testGeomData.integrationElements(point) *
                          trialGeomData.integrationElements(point) *
                          quadWeights[point];

    outOfPlaceConjugateTransposeDimAndDofDimensions(testValues[transIndex],
                                                    tmpTest);
    outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
        trialValues[transIndex], products, tmpTrial);

    arma::Mat<BasisFunctionType> matTest(&tmpTest[0], testDofCount,
                                         pointCount * transDim,
                                         false /* don't copy */, true);
    arma::Mat<BasisFunctionType> matTrial(&tmpTrial[0], trialDofCount,
                                          pointCount * transDim,
                                          false /* don't copy */, true);

    // this removes the complex conjugate from matTrial
    addABt(matTest, matTrial, result);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void evaluateWithTensorQuadratureRuleImpl(
    const GeometricalData<typename ScalarTraits<ResultType>::RealType> &
        testGeomData,
    const GeometricalData<typename ScalarTraits<ResultType>::RealType> &
        trialGeomData,
    const CollectionOf3dArrays<BasisFunctionType> &testValues,
    const CollectionOf3dArrays<BasisFunctionType> &trialValues,
    const CollectionOf4dArrays<KernelType> &kernelValues,
    const std::vector<typename ScalarTraits<ResultType>::RealType> &
        testQuadWeights,
    const std::vector<typename ScalarTraits<ResultType>::RealType> &
        trialQuadWeights,
    arma::Mat<ResultType> &result) {
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
  typedef typename Coercion<BasisFunctionType, ResultType>::Type
  IntermediateType;

  // Evaluate constants and assert that array dimensions are correct
  const size_t transCount = testValues.size();
  assert(transCount >= 1);
  assert(trialValues.size() == transCount);
  assert(kernelValues.size() == 1 || kernelValues.size() == transCount);

  const size_t testDofCount = testValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(testValues[i].extent(1) == testDofCount);
  const size_t trialDofCount = trialValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(trialValues[i].extent(1) == trialDofCount);

  const size_t testPointCount = testQuadWeights.size();
  const size_t trialPointCount = trialQuadWeights.size();

  for (size_t i = 0; i < kernelValues.size(); ++i) {
    assert(kernelValues[i].extent(0) == 1); // kernel is assumed to be scalar
    assert(kernelValues[i].extent(1) == 1); // kernel is assumed to be scalar
    assert(kernelValues[i].extent(2) == testPointCount);
    assert(kernelValues[i].extent(3) == trialPointCount);
  }
  for (size_t i = 0; i < transCount; ++i)
    assert(testValues[i].extent(2) == testPointCount);
  for (size_t i = 0; i < transCount; ++i)
    assert(trialValues[i].extent(2) == trialPointCount);

  assert(result.n_rows == testDofCount);
  assert(result.n_cols == trialDofCount);

  // Initialize the result matrix
  result.fill(0);

  // Declare temporary memory areas
  std::vector<ResultType, tbb::scalable_allocator<ResultType>> tmpReordered,
      tmpIntermediate;

  for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
    // Multiply elements of all test- and trialValues arrays and
    // transpose their dim and dof dimensions
    const size_t transDim = testValues[transIndex].extent(0);
    assert(transDim > 0);
    assert(trialValues[transIndex].extent(0) == transDim);

    const size_t kernelIndex = kernelValues.size() == 1 ? 0 : transIndex;
    arma::Mat<KernelType> matKernel(
        const_cast<KernelType *>(kernelValues[kernelIndex].begin()),
        testPointCount, trialPointCount, false /* don't copy */, true);

    if (testDofCount >= trialDofCount) {
      // TmpIntermediate = Kernel * Trial
      outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
          trialValues[transIndex], trialGeomData, trialQuadWeights,
          tmpReordered);
      arma::Mat<ResultType> matTrial(&tmpReordered[0], trialDofCount * transDim,
                                     trialPointCount, false);
      tmpIntermediate.resize(matTrial.n_rows * matKernel.n_rows);

      arma::Mat<IntermediateType> matTmp(
          &tmpIntermediate[0], trialDofCount * transDim, testPointCount, false);
      matTmp = matTrial * matKernel.t();
      matTmp.set_size(trialDofCount, transDim * testPointCount);

      // Result += Test * Tmp
      outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
          testValues[transIndex], testGeomData, testQuadWeights, tmpReordered);
      arma::Mat<ResultType> matTest(&tmpReordered[0], testDofCount,
                                    transDim * testPointCount, false);
      result += matTest * matTmp.t();
    } else { // testDofCount < trialDofCount
      // TmpIntermediate = Test * Kernel
      outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
          testValues[transIndex], testGeomData, testQuadWeights, tmpReordered);
      arma::Mat<ResultType> matTest(&tmpReordered[0], testDofCount * transDim,
                                    testPointCount, false);
      tmpIntermediate.resize(matTest.n_rows * matKernel.n_cols);

      arma::Mat<IntermediateType> matTmp(
          &tmpIntermediate[0], testDofCount * transDim, trialPointCount, false);
      matTmp = matTest * matKernel;
      matTmp.set_size(testDofCount, transDim * trialPointCount);

      // Result += Tmp * Trial
      outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
          trialValues[transIndex], trialGeomData, trialQuadWeights,
          tmpReordered);
      arma::Mat<ResultType> matTrial(&tmpReordered[0], trialDofCount,
                                     transDim * trialPointCount, false);
      result += matTmp * matTrial.t();
    }
  }
}

} // namespace

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void TypicalTestScalarKernelTrialIntegralBase<
    BasisFunctionType, KernelType,
    ResultType>::addGeometricalDependencies(size_t &testGeomDeps,
                                            size_t &trialGeomDeps) const {
  testGeomDeps |= INTEGRATION_ELEMENTS;
  trialGeomDeps |= INTEGRATION_ELEMENTS;
}

template <typename BasisFunctionType_, typename ResultType_>
void TypicalTestScalarKernelTrialIntegral<BasisFunctionType_,
                                          BasisFunctionType_, ResultType_>::
    evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf3dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &quadWeights,
        arma::Mat<ResultType> &result) const {
  evaluateWithNontensorQuadratureRuleStandardImpl(
      testGeomData, trialGeomData, testValues, trialValues, kernelValues,
      quadWeights, result);
}

template <typename CoordinateType_>
void TypicalTestScalarKernelTrialIntegral<CoordinateType_,
                                          std::complex<CoordinateType_>,
                                          std::complex<CoordinateType_>>::
    evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf4dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        arma::Mat<ResultType> &result) const {
  evaluateWithTensorQuadratureRuleImpl(
      testGeomData, trialGeomData, testValues, trialValues, kernelValues,
      testQuadWeights, trialQuadWeights, result);
}

template <typename CoordinateType>
void TypicalTestScalarKernelTrialIntegral<CoordinateType,
                                          std::complex<CoordinateType>,
                                          std::complex<CoordinateType>>::
    evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf3dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &quadWeights,
        arma::Mat<ResultType> &result) const {
  // We assume the integrand has the structure
  // (i -- transformations)
  // sum_i (\vec test_i \cdot \vec trial_i) kernel
  // or sum_i (\vec test_i \cdot kernel_i \vec trial_i) kernel

  // Evaluate constants and assert that array dimensions are correct
  const size_t transCount = testValues.size();
  assert(transCount >= 1);
  assert(trialValues.size() == transCount);
  assert(kernelValues.size() == 1 || kernelValues.size() == transCount);

  const size_t testDofCount = testValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(testValues[i].extent(1) == testDofCount);
  const size_t trialDofCount = trialValues[0].extent(1);
  for (size_t i = 1; i < transCount; ++i)
    assert(trialValues[i].extent(1) == trialDofCount);
  const size_t pointCount = quadWeights.size();

  for (size_t i = 0; i < kernelValues.size(); ++i)
    assert(kernelValues[i].extent(2) == pointCount);
  for (size_t i = 0; i < transCount; ++i)
    assert(testValues[i].extent(2) == pointCount);
  for (size_t i = 0; i < transCount; ++i)
    assert(trialValues[i].extent(2) == pointCount);
  assert(result.n_rows == testDofCount);
  assert(result.n_cols == trialDofCount);

  // Allocate memory for temporary arrays
  std::vector<CoordinateType, tbb::scalable_allocator<CoordinateType>>
  productsReal(pointCount);
  std::vector<CoordinateType, tbb::scalable_allocator<CoordinateType>>
  productsImag(pointCount);
  std::vector<CoordinateType, tbb::scalable_allocator<CoordinateType>> tmpTest;
  std::vector<CoordinateType, tbb::scalable_allocator<CoordinateType>> tmpTrial;
  arma::Mat<CoordinateType> matResultReal(testDofCount, trialDofCount);
  arma::Mat<CoordinateType> matResultImag(testDofCount, trialDofCount);

  // Evaluate each term of the integral in term
  for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
    const size_t transDim = testValues[transIndex].extent(0);
    assert(trialValues[transIndex].extent(0) == transDim);

    if (transIndex == 0 || kernelValues.size() > 1)
      for (size_t point = 0; point < pointCount; ++point) {
        CoordinateType partialProduct =
            testGeomData.integrationElements(point) *
            trialGeomData.integrationElements(point) * quadWeights[point];
        productsReal[point] =
            realPart(kernelValues[transIndex](0, 0, point)) * partialProduct;
        productsImag[point] =
            imagPart(kernelValues[transIndex](0, 0, point)) * partialProduct;
      }

    outOfPlaceConjugateTransposeDimAndDofDimensions(testValues[transIndex],
                                                    tmpTest);
    arma::Mat<BasisFunctionType> matTest(&tmpTest[0], testDofCount,
                                         pointCount * transDim,
                                         false /* don't copy */, true);

    // Process the real part of the kernel
    outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
        trialValues[transIndex], productsReal, tmpTrial);
    {
      arma::Mat<BasisFunctionType> matTrial(&tmpTrial[0], trialDofCount,
                                            pointCount * transDim,
                                            false /* don't copy */, true);
      if (transIndex == 0)
        matResultReal = matTest * matTrial.t();
      else
        matResultReal += matTest * matTrial.t();
    }

    // Process the imaginary part of the kernel
    outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
        trialValues[transIndex], productsImag, tmpTrial);
    {
      arma::Mat<BasisFunctionType> matTrial(&tmpTrial[0], trialDofCount,
                                            pointCount * transDim,
                                            false /* don't copy */, true);
      if (transIndex == 0)
        matResultImag = matTest * matTrial.t();
      else
        matResultImag += matTest * matTrial.t();
    }
  }
  result = arma::Mat<ResultType>(matResultReal, matResultImag);
}

template <typename BasisFunctionType_, typename ResultType_>
void TypicalTestScalarKernelTrialIntegral<BasisFunctionType_,
                                          BasisFunctionType_, ResultType_>::
    evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf4dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        arma::Mat<ResultType> &result) const {
  evaluateWithTensorQuadratureRuleImpl(
      testGeomData, trialGeomData, testValues, trialValues, kernelValues,
      testQuadWeights, trialQuadWeights, result);
}

template <typename CoordinateType>
TypicalTestScalarKernelTrialIntegral<
    std::complex<CoordinateType>, CoordinateType,
    std::complex<CoordinateType>>::TypicalTestScalarKernelTrialIntegral()
    : m_standardIntegral(TestScalarKernelTrialIntegrandFunctor<
          BasisFunctionType, KernelType, ResultType>()) {}

template <typename CoordinateType>
void TypicalTestScalarKernelTrialIntegral<std::complex<CoordinateType>,
                                          CoordinateType,
                                          std::complex<CoordinateType>>::
    evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf3dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &quadWeights,
        arma::Mat<ResultType> &result) const {
  evaluateWithNontensorQuadratureRuleStandardImpl(
      testGeomData, trialGeomData, testValues, trialValues, kernelValues,
      quadWeights, result);
}

template <typename CoordinateType>
void TypicalTestScalarKernelTrialIntegral<std::complex<CoordinateType>,
                                          CoordinateType,
                                          std::complex<CoordinateType>>::
    evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf4dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        arma::Mat<ResultType> &result) const {
  m_standardIntegral.evaluateWithTensorQuadratureRule(
      testGeomData, trialGeomData, testValues, trialValues, kernelValues,
      testQuadWeights, trialQuadWeights, result);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    TypicalTestScalarKernelTrialIntegral);

} // namespace Fiber

#endif
