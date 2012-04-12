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

#include "adjoint_double_layer_potential_3d_kernel.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>

#include "geometrical_data.hpp"

namespace Fiber
{

// Double potential: derivative wrt. trial normal

template <typename ValueType>
void AdjointDoubleLayerPotential3DKernel<ValueType>::addGeometricalDependencies(
        int& testGeomDeps, int& trialGeomDeps) const
{
    testGeomDeps |= GLOBALS | NORMALS;
    trialGeomDeps |= GLOBALS;
}

template <typename ValueType>
inline ValueType AdjointDoubleLayerPotential3DKernel<ValueType>::evaluateAtPointPair(
        const arma::Col<ValueType>& testPoint,
        const arma::Col<ValueType>& trialPoint,
        const arma::Col<ValueType>& testNormal) const
{
    const int coordCount = testPoint.n_rows;

    ValueType numeratorSum = 0., denominatorSum = 0.;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
    {
        ValueType diff = testPoint(coordIndex) - trialPoint(coordIndex);
        denominatorSum += diff * diff;
        numeratorSum += diff * testNormal(coordIndex);
    }
    ValueType distance = sqrt(denominatorSum);
    return -numeratorSum / (4. * M_PI * distance * distance * distance);
}

template <typename ValueType>
void AdjointDoubleLayerPotential3DKernel<ValueType>::evaluateAtPointPairs(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        arma::Cube<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;
    const arma::Mat<ValueType>& testNormals = testGeomData.normals;

#ifndef NDEBUG
    const int worldDim = worldDimension();
    if (testPoints.n_rows != worldDim || trialPoints.n_rows != worldDim)
        throw std::invalid_argument("AdjointDoubleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "3D coordinates required");
    if (testPoints.n_cols != trialPoints.n_cols)
        throw std::invalid_argument("AdjointDoubleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "number of test and trial points must be equal");
    assert(testNormals.n_rows == worldDim);
    assert(testNormals.n_cols == testPoints.n_cols);
#endif

    const int pointCount = testPoints.n_cols;
    result.set_size(1, 1, pointCount);
    for (int i = 0; i < pointCount; ++i)
        result(0, 0, i) = evaluateAtPointPair(
                    testPoints.unsafe_col(i), trialPoints.unsafe_col(i),
                    testNormals.unsafe_col(i));
}

template <typename ValueType>
void AdjointDoubleLayerPotential3DKernel<ValueType>::evaluateOnGrid(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        Array4D<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;
    const arma::Mat<ValueType>& testNormals = testGeomData.normals;

#ifndef NDEBUG
    const int worldDim = worldDimension();
    if (testPoints.n_rows != worldDim || trialPoints.n_rows != worldDim)
        throw std::invalid_argument("AdjointDoubleLayerPotential3DKernel::evaluate(): "
                                    "3D coordinates required");
    assert(testNormals.n_rows == worldDim);
    assert(testNormals.n_cols == testPoints.n_cols);
#endif

    const int testPointCount = testPoints.n_cols;
    const int trialPointCount = trialPoints.n_cols;
    result.set_size(1, testPointCount, 1, trialPointCount);
    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
            result(0, testIndex, 0, trialIndex) = evaluateAtPointPair(
                        testPoints.unsafe_col(testIndex),
                        trialPoints.unsafe_col(trialIndex),
                        testNormals.unsafe_col(testIndex));
}

template<typename ValueType>
std::pair<const char*,int> AdjointDoubleLayerPotential3DKernel<ValueType>::evaluateClCode () const
{
    // TODO!!!
    //return std::string (adjoint_double_layer_potential_3D_kernel_cl,
    //		adjoint_double_layer_potential_3D_kernel_cl_len);
    throw std::runtime_error ("AdjointDoubleLayerPotential3DKernel::evaluateClCode not implemented\n");
    return std::pair<const char*,int>();
}

#ifdef COMPILE_FOR_FLOAT
template class AdjointDoubleLayerPotential3DKernel<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class AdjointDoubleLayerPotential3DKernel<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class AdjointDoubleLayerPotential3DKernel<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class AdjointDoubleLayerPotential3DKernel<std::complex<double> >;
#endif

} // namespace Fiber
