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

#include "default_test_single_scalar_kernel_trial_integral.hpp"

#include "collection_of_4d_arrays.hpp"
#include "explicit_instantiation.hpp"
#include "geometrical_data.hpp"
#include "../common/acc.hpp"
#include "../common/complex_aux.hpp"

#include <cassert>
#include <iostream>

namespace Fiber
{

namespace
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void evaluateWithNontensorQuadratureRuleStandardImpl(
        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& testGeomData,
        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<typename ScalarTraits<ResultType>::RealType>& quadWeights,
        arma::Mat<ResultType>& result)
{
    // std::cout << "var 1" << std::endl;
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

    std::vector<KernelType> products(pointCount);
    std::vector<BasisFunctionType> tmp;

    for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
        const size_t transDim = testValues[transIndex].extent(0);
        assert(trialValues[transIndex].extent(0) == transDim);

        if (transIndex == 0 || kernelValues.size() > 1)
            for (size_t point = 0; point < pointCount; ++point)
                products[point] = kernelValues[0](0, 0, point) *
                        testGeomData.integrationElements(point) *
                        trialGeomData.integrationElements(point) *
                        quadWeights[point];

        tmp.resize(testDofCount * transDim);
        // testValues[transIndex].set_size(testDofCount, transDim, pointCount);
        for (size_t point = 0; point < pointCount; ++point) {
            memcpy(&tmp[0],
                   testValues[transIndex].begin() + point * tmp.size(),
                   tmp.size() * sizeof(BasisFunctionType));
            arma::Mat<BasisFunctionType> origMat(
                        &tmp[0], transDim, testDofCount, false, true);
            arma::Mat<BasisFunctionType> transposedMat(
                        testValues[transIndex].begin() + point * tmp.size(),
                        testDofCount, transDim, false, true);
            transposedMat = origMat.t(); // we take the complex conjugate here
        }

        tmp.resize(trialDofCount * transDim);
        // trialValues[transIndex].set_size(trialDofCount, transDim, pointCount);
        for (size_t point = 0; point < pointCount; ++point) {
            memcpy(&tmp[0],
                   trialValues[transIndex].begin() + point * tmp.size(),
                   tmp.size() * sizeof(BasisFunctionType));
            arma::Mat<BasisFunctionType> origMat(
                        &tmp[0], transDim, trialDofCount, false, true);
            arma::Mat<BasisFunctionType> transposedMat(
                        trialValues[transIndex].begin() + point * tmp.size(),
                        trialDofCount, transDim, false, true);
            // we take the complex conj. here, afterwards we'll remove it
            transposedMat = conj(products[point]) * origMat.t();
        }

        arma::Mat<BasisFunctionType> matTest(testValues[transIndex].begin(),
                                             testDofCount, pointCount  * transDim,
                                             false /* don't copy */, true);
        arma::Mat<BasisFunctionType> matTrial(trialValues[transIndex].begin(),
                                              trialDofCount, pointCount * transDim,
                                              false /* don't copy */, true);

        // this removes the complex conjugate from matTrial
        result = result + matTest * matTrial.t();
    }
}

template <typename BasisFunctionType>
void inplaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
        const GeometricalData<typename ScalarTraits<BasisFunctionType>::RealType>& geomData,
        const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>& quadWeights,
        std::vector<BasisFunctionType>& tmp,
        _3dArray<BasisFunctionType>& values)
{
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    const size_t transDim = values.extent(0);
    const size_t dofCount = values.extent(1);
    const size_t pointCount = values.extent(2);
    const size_t sliceSize = transDim * dofCount;
    if (transDim > 1)
        tmp.resize(dofCount * transDim);
    // First, transpose logical dimensions...
    //values.set_size(dofCount, transDim, pointCount);
    BasisFunctionType* sliceStart = values.begin();
    for (size_t point = 0; point < pointCount; ++point) {
        arma::Mat<BasisFunctionType> transposedMat(
                    sliceStart, dofCount, transDim, false, true);
        const CoordinateType weight =
                geomData.integrationElements(point) * quadWeights[point];
        // And now, transpose data physically (if needed, i.e. if slices
        // are not logically one-dimensional)
        if (transDim > 1) {
            memcpy(&tmp[0], sliceStart, sliceSize * sizeof(BasisFunctionType));
            arma::Mat<BasisFunctionType> origMat(
                        &tmp[0], transDim, dofCount, false, true);
            // we take the complex conjugate here
            transposedMat = weight * origMat.t();
        } else // no need to transpose data physically
            transposedMat *= weight;
        sliceStart += sliceSize;
    }
}

template <typename BasisFunctionType, typename ResultType>
void setToProduct(const arma::Mat<BasisFunctionType>& source,
                  typename ScalarTraits<BasisFunctionType>::RealType weight,
                  arma::Mat<ResultType>& result)
{
    result = weight * source.t();
}

template <typename CoordinateType>
void setToProduct(const arma::Mat<CoordinateType>& source,
                  CoordinateType weight,
                  arma::Mat<std::complex<CoordinateType> >& result)
{
    result.set_real(weight * source.t());
}

template <typename BasisFunctionType, typename ResultType>
void outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
        const _3dArray<BasisFunctionType>& values,
        const GeometricalData<typename ScalarTraits<BasisFunctionType>::RealType>& geomData,
        const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>& quadWeights,
        std::vector<ResultType>& tmp)
{
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    const size_t transDim = values.extent(0);
    const size_t dofCount = values.extent(1);
    const size_t pointCount = values.extent(2);
    const size_t sliceSize = transDim * dofCount;

    const BasisFunctionType* sourceSliceStart = values.begin();
    tmp.resize(dofCount * transDim * pointCount);
    ResultType* destSliceStart = &tmp[0];

    for (size_t point = 0; point < pointCount; ++point) {
        arma::Mat<BasisFunctionType> origMat(
                    const_cast<BasisFunctionType*>(sourceSliceStart),
                    transDim, dofCount, false);
        arma::Mat<ResultType> transposedMat(
                    destSliceStart,
                    dofCount, transDim, false);
        const CoordinateType weight =
                geomData.integrationElements(point) * quadWeights[point];
        // And now, transpose data physically (if needed, i.e. if slices
        // are not logically one-dimensional)
        // we take the complex conjugate here
        setToProduct(origMat, weight, transposedMat);
        // transposedMat = weight * origMat.t();
        sourceSliceStart += sliceSize;
        destSliceStart += sliceSize;
    }
}

template <typename BasisFunctionType, typename ResultType>
inline void addABt(
        const arma::Mat<BasisFunctionType>& A,
        const arma::Mat<BasisFunctionType>& B,
        arma::Mat<ResultType>& C)
{
    C = C + A * B.t();
}

template <typename ResultType>
inline void addABt(
        const arma::Mat<ResultType>& A,
        const arma::Mat<ResultType>& B,
        arma::Mat<ResultType>& C)
{
    C += A * B.t();
}

template <typename BasisFunctionType, typename ResultType>
inline void multiplyTestKernelTrialAndAddToResult(
        const arma::Mat<BasisFunctionType>& matTest,
        const arma::Mat<BasisFunctionType>& matKernel,
        const arma::Mat<BasisFunctionType>& matTrial,
        std::vector<BasisFunctionType>& tmp,
        arma::Mat<ResultType>& result)
{
    tmp.resize(matTrial.n_rows * matKernel.n_rows);
    arma::Mat<BasisFunctionType> matTmp(&tmp[0], matTrial.n_rows, matKernel.n_rows,
                                 false, true);
    matTmp = matTrial * matKernel.t();
    addABt(matTest, matTmp, result);
}

//template <typename BasisFunctionType, typename ResultType>
//void evaluateWithTensorQuadratureRuleImpl(
//        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& testGeomData,
//        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& trialGeomData,
//        CollectionOf3dArrays<BasisFunctionType>& testValues,
//        CollectionOf3dArrays<BasisFunctionType>& trialValues,
//        CollectionOf4dArrays<BasisFunctionType>& kernelValues,
//        const std::vector<typename ScalarTraits<ResultType>::RealType>& testQuadWeights,
//        const std::vector<typename ScalarTraits<ResultType>::RealType>& trialQuadWeights,
//        arma::Mat<ResultType>& result)
//{
//    typedef BasisFunctionType KernelType;
//    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

//    // Evaluate constants and assert that array dimensions are correct
//    const size_t transCount = testValues.size();
//    assert(transCount >= 1);
//    assert(trialValues.size() == transCount);
//    assert(kernelValues.size() == 1 || kernelValues.size() == transCount);

//    const size_t testDofCount = testValues[0].extent(1);
//    for (size_t i = 1; i < transCount; ++i)
//        assert(testValues[i].extent(1) == testDofCount);
//    const size_t trialDofCount = trialValues[0].extent(1);
//    for (size_t i = 1; i < transCount; ++i)
//        assert(trialValues[i].extent(1) == trialDofCount);

//    const size_t testPointCount = testQuadWeights.size();
//    const size_t trialPointCount = trialQuadWeights.size();

//    for (size_t i = 0; i < kernelValues.size(); ++i) {
//        assert(kernelValues[i].extent(0) == 1); // kernel is assumed to be scalar
//        assert(kernelValues[i].extent(1) == 1); // kernel is assumed to be scalar
//        assert(kernelValues[i].extent(2) == testPointCount);
//        assert(kernelValues[i].extent(3) == trialPointCount);
//    }
//    for (size_t i = 0; i < transCount; ++i)
//        assert(testValues[i].extent(2) == testPointCount);
//    for (size_t i = 0; i < transCount; ++i)
//        assert(trialValues[i].extent(2) == trialPointCount);

//    assert(result.n_rows == testDofCount);
//    assert(result.n_cols == trialDofCount);

//    // Initialize the result matrix
//    result.fill(0);

//    // Declare temporary memory area
//    std::vector<BasisFunctionType> tmp;

//    for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
//        // Multiply elements of all test- and trialValues arrays and
//        // transpose their dim and dof dimensions
//        const size_t transDim = testValues[transIndex].extent(0);
//        assert(transDim > 0);
//        assert(trialValues[transIndex].extent(0) == transDim);
//        inplaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
//                    testGeomData, testQuadWeights, tmp, testValues[transIndex]);
//        inplaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
//                    trialGeomData, trialQuadWeights, tmp, trialValues[transIndex]);

//        const size_t kernelIndex = kernelValues.size() == 1 ? 0 : transIndex;
//        // Now do the matrix-matrix-matrix multiplications
//        arma::Mat<BasisFunctionType> matTest(
//                    testValues[transIndex].begin(),
//                    testDofCount, testPointCount * transDim,
//                    false /* don't copy */, true);
//        arma::Mat<KernelType> matKernel(
//                    kernelValues[kernelIndex].begin(),
//                    testPointCount, trialPointCount,
//                    false /* don't copy */, true);
//        arma::Mat<BasisFunctionType> matTrial(
//                    trialValues[transIndex].begin(),
//                    trialDofCount, trialPointCount * transDim,
//                    false /* don't copy */, true);
//        multiplyTestKernelTrialAndAddToResult(
//                    matTest, matKernel, matTrial, tmp, result);
//    }
//}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void evaluateWithTensorQuadratureRuleImpl(
        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& testGeomData,
        const GeometricalData<typename ScalarTraits<ResultType>::RealType>& trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType>& testValues,
        const CollectionOf3dArrays<BasisFunctionType>& trialValues,
        const CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<typename ScalarTraits<ResultType>::RealType>& testQuadWeights,
        const std::vector<typename ScalarTraits<ResultType>::RealType>& trialQuadWeights,
        arma::Mat<ResultType>& result)
{
    // typedef ResultType KernelType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    typedef typename Coercion<BasisFunctionType, ResultType>::Type IntermediateType;

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
    std::vector<ResultType> tmpReordered, tmpIntermediate;

    for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
        // Multiply elements of all test- and trialValues arrays and
        // transpose their dim and dof dimensions
        const size_t transDim = testValues[transIndex].extent(0);
        assert(transDim > 0);
        assert(trialValues[transIndex].extent(0) == transDim);

        const size_t kernelIndex = kernelValues.size() == 1 ? 0 : transIndex;
        arma::Mat<KernelType> matKernel(
                    const_cast<KernelType*>(kernelValues[kernelIndex].begin()),
                    testPointCount, trialPointCount,
                    false /* don't copy */, true);

        if (testDofCount >= trialDofCount) {
            // TmpIntermediate = Kernel * Trial
            outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
                        trialValues[transIndex], trialGeomData, trialQuadWeights,
                        tmpReordered);
            arma::Mat<ResultType> matTrial(
                        &tmpReordered[0],
                        trialDofCount * transDim, trialPointCount, false);
            tmpIntermediate.resize(matTrial.n_rows * matKernel.n_rows);

            arma::Mat<IntermediateType> matTmp(
                        &tmpIntermediate[0],
                        trialDofCount * transDim, testPointCount, false);
            matTmp = matTrial * matKernel.t();
            matTmp.set_size(trialDofCount, transDim * testPointCount);

            // Result += Test * Tmp
            outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
                        testValues[transIndex], testGeomData, testQuadWeights,
                        tmpReordered);
            arma::Mat<ResultType> matTest(
                        &tmpReordered[0],
                        testDofCount, transDim * testPointCount, false);
            result += matTest * matTmp.t();
        } else { // testDofCount < trialDofCount
            // TmpIntermediate = Test * Kernel
            outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
                        testValues[transIndex], testGeomData, testQuadWeights,
                        tmpReordered);
            arma::Mat<ResultType> matTest(
                        &tmpReordered[0],
                        testDofCount * transDim, testPointCount, false);
            tmpIntermediate.resize(matTest.n_rows * matKernel.n_cols);

            arma::Mat<IntermediateType> matTmp(
                        &tmpIntermediate[0],
                        testDofCount * transDim, trialPointCount, false);
            matTmp = matTest * matKernel;
            matTmp.set_size(testDofCount, transDim * trialPointCount);

            // Result += Tmp * Trial
            outOfPlaceMultiplyByWeightsAndConjugateTransposeDimAndDofDimensions(
                        trialValues[transIndex], trialGeomData, trialQuadWeights,
                        tmpReordered);
            arma::Mat<ResultType> matTrial(
                        &tmpReordered[0],
                        trialDofCount, transDim * trialPointCount, false);
            result += matTmp * matTrial.t();
        }
    }
}

} // namespace

template <typename BasisFunctionType_, typename ResultType_>
void DefaultTestSingleScalarKernelTrialIntegral<BasisFunctionType_,
BasisFunctionType_, ResultType_>::
evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& quadWeights,
        arma::Mat<ResultType>& result) const
{
    evaluateWithNontensorQuadratureRuleStandardImpl(
                testGeomData, trialGeomData,
                testValues, trialValues, kernelValues, quadWeights, result);
}

template <typename CoordinateType>
void DefaultTestSingleScalarKernelTrialIntegral<CoordinateType,
std::complex<CoordinateType>, std::complex<CoordinateType> >::
evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& quadWeights,
        arma::Mat<ResultType>& result) const
{
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
    std::vector<KernelType> products(pointCount);
    arma::Mat<CoordinateType> matResultReal(testDofCount, trialDofCount);
    arma::Mat<CoordinateType> matResultImag(testDofCount, trialDofCount);

    // Evaluate each term of the integral in term
    for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
        const size_t transDim = testValues[transIndex].extent(0);
        assert(trialValues[transIndex].extent(0) == transDim);

        if (transIndex == 0 || kernelValues.size() > 1)
            for (size_t point = 0; point < pointCount; ++point)
                products[point] = kernelValues[0](0, 0, point) *
                        testGeomData.integrationElements(point) *
                        trialGeomData.integrationElements(point) *
                        quadWeights[point];

        // this loop can be removed for transDim == 1, because the internal
        // transpose is not needed;
        std::vector<CoordinateType> tmp;
        tmp.resize(testDofCount * transDim);
        // transpose first two dimensions: first logically...
        //testValues[transIndex].set_size(testDofCount, transDim, pointCount);
        BasisFunctionType* sliceStart = testValues[transIndex].begin();
        for (size_t point = 0; point < pointCount; ++point) {
            memcpy(&tmp[0], sliceStart, tmp.size() * sizeof(BasisFunctionType));
            arma::Mat<BasisFunctionType> source(
                        &tmp[0], transDim, testDofCount, false, true);
            arma::Mat<BasisFunctionType> dest(
                        sliceStart, testDofCount, transDim, false, true);
            // and then physically:
            dest = source.t();
            sliceStart += testDofCount * transDim;
        }
        arma::Mat<BasisFunctionType> matTest(testValues[transIndex].begin(),
                                             testDofCount, pointCount  * transDim,
                                             false /* don't copy */, true);

        // Process the real part of the kernel
        std::vector<CoordinateType> scaledTrialValues;
        scaledTrialValues.resize(trialDofCount * transDim * pointCount);
        {
            BasisFunctionType* sourceSliceStart =
                    trialValues[transIndex].begin();
            BasisFunctionType* destSliceStart = &scaledTrialValues[0];
            for (size_t point = 0; point < pointCount; ++point) {
                arma::Mat<BasisFunctionType> source(
                            sourceSliceStart,
                            transDim, trialDofCount, false, true);
                arma::Mat<BasisFunctionType> dest(
                            destSliceStart,
                            trialDofCount, transDim, false, true);
//                dest = realPart(products[point]) * source.t();
                dest = source.t();
                dest *= realPart(products[point]);
                sourceSliceStart += trialDofCount * transDim;
                destSliceStart += trialDofCount * transDim;
            }

            arma::Mat<BasisFunctionType> matTrial(&scaledTrialValues[0],
                                                  trialDofCount, pointCount * transDim,
                                                  false /* don't copy */, true);
            if (transIndex == 0)
                matResultReal = matTest * matTrial.t();
            else
                matResultReal += matTest * matTrial.t();
        }

        // Now process the imaginary part of the kernel
        {
            BasisFunctionType* sourceSliceStart =
                    trialValues[transIndex].begin();
            BasisFunctionType* destSliceStart = &scaledTrialValues[0];
            for (size_t point = 0; point < pointCount; ++point) {
                arma::Mat<BasisFunctionType> source(
                            sourceSliceStart,
                            transDim, trialDofCount, false, true);
                arma::Mat<BasisFunctionType> dest(
                            destSliceStart,
                            trialDofCount, transDim, false, true);
                dest = source.t();
                dest *= imagPart(products[point]);
                sourceSliceStart += trialDofCount * transDim;
                destSliceStart += trialDofCount * transDim;
            }

            arma::Mat<BasisFunctionType> matTrial(&scaledTrialValues[0],
                                                  trialDofCount, pointCount * transDim,
                                                  false /* don't copy */, true);
            if (transIndex == 0)
                matResultImag = matTest * matTrial.t();
            else
                matResultImag += matTest * matTrial.t();
        }
    }
    result = arma::Mat<ResultType>(matResultReal, matResultImag);
}

template <typename CoordinateType>
void DefaultTestSingleScalarKernelTrialIntegral<std::complex<CoordinateType>,
CoordinateType, std::complex<CoordinateType> >::
evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& quadWeights,
        arma::Mat<ResultType>& result) const
{
    evaluateWithNontensorQuadratureRuleStandardImpl(
                testGeomData, trialGeomData,
                testValues, trialValues, kernelValues, quadWeights, result);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void DefaultTestSingleScalarKernelTrialIntegralBase<BasisFunctionType, KernelType, ResultType>::
addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const
{
    testGeomDeps |= INTEGRATION_ELEMENTS;
    trialGeomDeps |= INTEGRATION_ELEMENTS;

    m_functor.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
}

//template <typename BasisFunctionType, typename KernelType, typename ResultType>
//void DefaultTestSingleScalarKernelTrialIntegralBase<BasisFunctionType, KernelType, ResultType>::
//evaluateWithTensorQuadratureRule(
//        const GeometricalData<CoordinateType>& testGeomData,
//        const GeometricalData<CoordinateType>& trialGeomData,
//        const CollectionOf3dArrays<BasisFunctionType>& testValues,
//        const CollectionOf3dArrays<BasisFunctionType>& trialValues,
//        const CollectionOf4dArrays<KernelType>& kernelValues,
//        const std::vector<CoordinateType>& testQuadWeights,
//        const std::vector<CoordinateType>& trialQuadWeights,
//        arma::Mat<ResultType>& result) const
//{
//    // Evaluate constants

//    const size_t testDofCount = testValues[0].extent(1);
//    const size_t trialDofCount = trialValues[0].extent(1);

//    const size_t testPointCount = testQuadWeights.size();
//    const size_t trialPointCount = trialQuadWeights.size();

//    // Assert that array dimensions are correct

//    for (size_t i = 0; i < kernelValues.size(); ++i) {
//        assert(kernelValues[i].extent(2) == testPointCount);
//        assert(kernelValues[i].extent(3) == trialPointCount);
//    }
//    for (size_t i = 0; i < testValues.size(); ++i)
//        assert(testValues[i].extent(2) == testPointCount);
//    for (size_t i = 0; i < trialValues.size(); ++i)
//        assert(trialValues[i].extent(2) == trialPointCount);

//    // Integrate

//    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
//        for (size_t testDof = 0; testDof < testDofCount; ++testDof)
//        {
//            ResultType sum = 0.;
//            for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
//                const CoordinateType trialWeight =
//                        trialGeomData.integrationElements(trialPoint) *
//                        trialQuadWeights[trialPoint];
//                ResultType partialSum = 0.;
//                for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
//                    const CoordinateType testWeight =
//                            testGeomData.integrationElements(testPoint) *
//                            testQuadWeights[testPoint];
//                    partialSum += m_functor.evaluate(
//                                testGeomData.const_slice(testPoint),
//                                trialGeomData.const_slice(trialPoint),
//                                testValues.const_slice(testDof, testPoint),
//                                trialValues.const_slice(trialDof, trialPoint),
//                                kernelValues.const_slice(testPoint, trialPoint)) *
//                            testWeight;
//                }
//                sum += partialSum * trialWeight;
//            }
//            result(testDof, trialDof) = sum;
//        }
//}

template <typename CoordinateType_>
void DefaultTestSingleScalarKernelTrialIntegral<CoordinateType_,
std::complex<CoordinateType_>, std::complex<CoordinateType_> >::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        arma::Mat<ResultType>& result) const
{
    evaluateWithTensorQuadratureRuleImpl(
                testGeomData, trialGeomData,
                testValues, trialValues, kernelValues,
                testQuadWeights, trialQuadWeights,
                result);
}

template <typename BasisFunctionType_, typename ResultType_>
void DefaultTestSingleScalarKernelTrialIntegral<BasisFunctionType_,
BasisFunctionType_, ResultType_>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        arma::Mat<ResultType>& result) const
{
    evaluateWithTensorQuadratureRuleImpl(
                testGeomData, trialGeomData,
                testValues, trialValues, kernelValues,
                testQuadWeights, trialQuadWeights,
                result);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void DefaultTestSingleScalarKernelTrialIntegralBase<BasisFunctionType, KernelType, ResultType>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        arma::Mat<ResultType>& result) const
{
    // throw std::runtime_error("Old implementation");
    // Evaluate constants and assert that array dimensions are correct
    const size_t transCount = testValues.size();
    assert(trialValues.size() == transCount);
    assert(kernelValues.size() == 1 ||
           kernelValues.size() == transCount);

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

    // Integrate

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (size_t testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
                const CoordinateType trialWeight =
                        trialGeomData.integrationElements(trialPoint) *
                        trialQuadWeights[trialPoint];
                ResultType partialSum = 0.;
                for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
                    const CoordinateType testWeight =
                            testGeomData.integrationElements(testPoint) *
                            testQuadWeights[testPoint];
                    partialSum += m_functor.evaluate(
                                testGeomData.const_slice(testPoint),
                                trialGeomData.const_slice(trialPoint),
                                testValues.const_slice(testDof, testPoint),
                                trialValues.const_slice(trialDof, trialPoint),
                                kernelValues.const_slice(testPoint, trialPoint)) *
                            testWeight;
                }
                sum += partialSum * trialWeight;
            }
            result(testDof, trialDof) = sum;
        }
}

//template <typename IntegrandFunctor>
//void DefaultTestSingleScalarKernelTrialIntegral<IntegrandFunctor>::
//evaluateWithNontensorQuadratureRule(
//        GeometricalData<CoordinateType>& testGeomData,
//        GeometricalData<CoordinateType>& trialGeomData,
//        CollectionOf3dArrays<BasisFunctionType>& testValues,
//        CollectionOf3dArrays<BasisFunctionType>& trialValues,
//        CollectionOf3dArrays<KernelType>& kernelValues,
//        std::vector<CoordinateType>& quadWeights,
//        arma::Mat<ResultType>& result) const
//{
//    // We assume the integrand has the structure
//    // (\vec test \cdot \vec trial_i) kernel
//    assert(testValues.size() == 1);
//    assert(trialValues.size() == 1);

//    // Evaluate constants
//    const size_t testDofCount = testValues[0].extent(1);
//    const size_t trialDofCount = trialValues[0].extent(1);
//    const size_t transDim = testValues[0].extent(0);
//    assert(trialValues[0].extent(0) == transDim);
//    const size_t pointCount = quadWeights.size();

////    if (transDim != 1)
////        throw std::runtime_error("transDim != 1");

//    // Assert that array dimensions are correct

//    assert(kernelValues[0].extent(2) == pointCount);
//    assert(testValues[0].extent(2) == pointCount);
//    assert(trialValues[0].extent(2) == pointCount);
//    assert(result.n_rows == testDofCount);
//    assert(result.n_cols == trialDofCount);

//    // Integrate

////    ProductDataReference products = m_products.local();
////    products.resize(pointCount);
////    for (size_t point = 0; point < pointCount; ++point)
////        products[point] = kernelValues[0](0, 0, point) *
////                testGeomData.integrationElements(point) *
////                trialGeomData.integrationElements(point) *
////                quadWeights[point];

////    FactorDataReference reorderedTestValues = m_reorderedTestValues.local();
////    reorderedTestValues.resize(transDim * testDofCount * pointCount);
////    for (size_t f = 0, r = 0; f < testDofCount; ++f)
////        for (size_t d = 0; d < transDim; ++d)
////            for (size_t p = 0; p < pointCount; ++p, ++r)
////                // Bempp::acc(reorderedTestValues, r) = testValues[0](d, f, p);
////                reorderedTestValues[r] = testValues[0](d, f, p);

////    FactorDataReference scaledTrialValues = m_scaledTrialValues.local();
////    scaledTrialValues.resize(transDim * trialDofCount * pointCount);
////    for (size_t f = 0, r = 0; f < trialDofCount; ++f)
////        for (size_t d = 0; d < transDim; ++d)
////            for (size_t p = 0; p < pointCount; ++p, ++r)
////                // Bempp::acc(scaledTrialValues, r) = trialValues[0](d, f, p) * products[p];
////                scaledTrialValues[r] = trialValues[0](d, f, p) * products[p];

////    arma::Mat<ResultType> matTest(&reorderedTestValues[0],
////                                   pointCount * transDim, testDofCount,
////                                   false /* don't copy */);
////    arma::Mat<ResultType> matTrial(&scaledTrialValues[0],
////                                   pointCount * transDim, trialDofCount,
////                                   false /* don't copy */);
////    // Teuchos::BLAS<int, ResultType> blas;
////    // blas.GEMM(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
////    //           testDofCount, trialDofCount,
////    //           pointCount * transDim,
////    //           1., &reorderedTestValues[0], 1, &scaledTrialValues[0], 1,
////    //           0., &result[0], 1.);

////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;
////    // std::cout << matTest.n_rows << " " << matTest.n_cols << std::endl;
////    // std::cout << matTrial.n_rows << " " << matTrial.n_cols << std::endl;
////    result = matTest.t() * matTrial;
////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;




////    ProductDataReference products = m_products.local();
////    products.resize(pointCount);
////    for (size_t point = 0; point < pointCount; ++point)
////        products[point] = kernelValues[0](0, 0, point) *
////                testGeomData.integrationElements(point) *
////                trialGeomData.integrationElements(point) *
////                quadWeights[point];

////    FactorDataReference reorderedTestValues = m_reorderedTestValues.local();
////    reorderedTestValues.resize(transDim * testDofCount * pointCount);
////    for (size_t p = 0, r = 0; p < pointCount; ++p)
////        for (size_t f = 0; f < testDofCount; ++f)
////            for (size_t d = 0; d < transDim; ++d, ++r)
////                // Bempp::acc(scaledTrialValues, r) = trialValues[0](d, f, p) * products[p];
////                reorderedTestValues[r] = testValues[0](d, f, p);

//////    FactorDataReference scaledTrialValues = m_scaledTrialValues.local();
//////    scaledTrialValues.resize(transDim * trialDofCount * pointCount);
//////    for (size_t p = 0, r = 0; p < pointCount; ++p)
//////        for (size_t f = 0; f < trialDofCount; ++f)
//////            for (size_t d = 0; d < transDim; ++d, ++r)
//////                // Bempp::acc(scaledTrialValues, r) = trialValues[0](d, f, p) * products[p];
//////                scaledTrialValues[r] = trialValues[0](d, f, p) * products[p];
////    FactorDataReference scaledTrialValues = m_scaledTrialValues.local();
////    scaledTrialValues.resize(transDim * trialDofCount * pointCount);
//////    for (size_t p = 0, r = 0; p < pointCount; ++p)
//////        for (size_t f = 0; f < trialDofCount; ++f)
//////            for (size_t d = 0; d < transDim; ++d, ++r)
//////                // Bempp::acc(scaledTrialValues, r) = trialValues[0](d, f, p) * products[p];
//////                scaledTrialValues[r] = trialValues[0](d, f, p) * products[p];
////    memcpy(&scaledTrialValues[0], trialValues[0].begin(), sizeof(CoordinateType) * scaledTrialValues.size());
////    for (size_t p = 0, r = 0; p < pointCount; ++p)
////        for (size_t f = 0; f < trialDofCount; ++f, ++r)
////                scaledTrialValues[r] *= products[p];

////    arma::Mat<ResultType> matTest(&reorderedTestValues[0],
////                                   testDofCount, pointCount,
////                                   false /* don't copy */);
////    arma::Mat<ResultType> matTrial(&scaledTrialValues[0],
////                                   trialDofCount, pointCount,
////                                   false /* don't copy */);

////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;
////    // std::cout << matTest.n_rows << " " << matTest.n_cols << std::endl;
////    // std::cout << matTrial.n_rows << " " << matTrial.n_cols << std::endl;
////    result = matTest * matTrial.t();// sctually conjugate would be needed...
////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;



////    ProductDataReference products = m_products.local();
////    products.resize(pointCount);
////    for (size_t point = 0; point < pointCount; ++point)
////        products[point] = kernelValues[0](0, 0, point) *
////                testGeomData.integrationElements(point) *
////                trialGeomData.integrationElements(point) *
////                quadWeights[point];
////    KernelType* pproducts = &products[0];

////    for (size_t p = 0, r = 0; p < pointCount; ++p)
////        for (size_t f = 0; f < trialDofCount; ++f)
////            for (size_t d = 0; d < transDim; ++d, ++r)
////                // Bempp::acc(scaledTrialValues, r) = trialValues[0](d, f, p) * products[p];
////                trialValues[0](d, f, p) *= pproducts[p];

//////    FactorDataReference scaledTrialValues = m_scaledTrialValues.local();
////    arma::Mat<BasisFunctionType> matTest(testValues[0].begin(),
////                                  testDofCount, pointCount,
////                                  false /* don't copy */);
////    arma::Mat<BasisFunctionType> matTrial(trialValues[0].begin(),
////                                   trialDofCount, pointCount,
////                                   false /* don't copy */);

////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;
////    // std::cout << matTest.n_rows << " " << matTest.n_cols << std::endl;
////    // std::cout << matTrial.n_rows << " " << matTrial.n_cols << std::endl;
////    arma::Mat<BasisFunctionType> resultX = matTest * matTrial.t();// sctually conjugate would be needed...
////    for (size_t i = 0; i < result.n_elem; ++i)
////        result[i] = resultX[i];
////    // std::cout << result.n_rows << " " << result.n_cols << std::endl;

////    arma::Mat<ResultType> compareResult(testDofCount, trialDofCount);
////    // std::cout << compareResult.n_rows << " " << compareResult.n_cols << std::endl;
////    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
////        for (size_t testDof = 0; testDof < testDofCount; ++testDof)
////        {
////            ResultType sum = 0.;
////            for (size_t point = 0; point < pointCount; ++point)
////                sum += m_functor.evaluate(
////                            testGeomData.const_slice(point),
////                            trialGeomData.const_slice(point),
////                            testValues.const_slice(testDof, point),
////                            trialValues.const_slice(trialDof, point),
////                            kernelValues.const_slice(point)) *
////                        (testGeomData.integrationElements(point) *
////                         trialGeomData.integrationElements(point) *
////                         quadWeights[point]);
//////            result(testDof, trialDof) = sum;
////            compareResult(testDof, trialDof) = sum;
////        }

//    ProductDataReference products = m_products.local();
//    products.resize(pointCount);
//    for (size_t point = 0; point < pointCount; ++point)
//        products[point] = kernelValues[0](0, 0, point) *
//                testGeomData.integrationElements(point) *
//                trialGeomData.integrationElements(point) *
//                quadWeights[point];

//    BFTDataReference bftData = m_bftData.local();
//    bftData.resize(testDofCount * transDim);
//    testValues[0].set_size(testDofCount, transDim, pointCount);
//    for (size_t point = 0; point < pointCount; ++point) {
//        memcpy(&bftData[0],
//               testValues[0].begin() + point * bftData.size(),
//               bftData.size() * sizeof(BasisFunctionType));
//        arma::Mat<BasisFunctionType> origMat(
//                    &bftData[0], transDim, testDofCount, false, true);
//        arma::Mat<BasisFunctionType> transposedMat(
//                    testValues[0].begin() + point * bftData.size(),
//                    testDofCount, transDim, false, true);
//        transposedMat = origMat.t();
//    }

//    bftData.resize(trialDofCount * transDim);
//    trialValues[0].set_size(trialDofCount, transDim, pointCount);
//    for (size_t point = 0; point < pointCount; ++point) {
//        memcpy(&bftData[0],
//               trialValues[0].begin() + point * bftData.size(),
//               bftData.size() * sizeof(BasisFunctionType));
//        arma::Mat<BasisFunctionType> origMat(
//                    &bftData[0], transDim, trialDofCount, false, true);
//        arma::Mat<BasisFunctionType> transposedMat(
//                    trialValues[0].begin() + point * bftData.size(),
//                    trialDofCount, transDim, false, true);
//        transposedMat = origMat.t();
//        transposedMat *= realPart(products[point]);
//    }

//    arma::Mat<BasisFunctionType> matTest(testValues[0].begin(),
//                                         testDofCount, pointCount  * transDim,
//                                  false /* don't copy */, true);
//    arma::Mat<BasisFunctionType> matTrial(trialValues[0].begin(),
//                                          trialDofCount, pointCount * transDim,
//                                          false /* don't copy */, true);

//    // std::cout << result.n_rows << " " << result.n_cols << std::endl;
//    // std::cout << matTest.n_rows << " " << matTest.n_cols << std::endl;
//    // std::cout << matTrial.n_rows << " " << matTrial.n_cols << std::endl;
//    arma::Mat<BasisFunctionType> resultX = matTest * matTrial.t();// sctually conjugate would be needed...
//    for (size_t i = 0; i < result.n_elem; ++i)
//        result[i] = resultX[i];

////    compareResult -= result;
////    std::cout << arma::norm(compareResult, 2) << std::endl;
//}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
        DefaultTestSingleScalarKernelTrialIntegral);

} // namespace Fiber

#endif
