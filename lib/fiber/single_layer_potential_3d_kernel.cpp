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

#include "single_layer_potential_3d_kernel.hpp"

#include <armadillo>
#include <cmath>

#include "geometrical_data.hpp"
#include "CL/single_layer_potential_3D_kernel.cl.str"

namespace Fiber
{

template <typename ValueType>
void SingleLayerPotential3DKernel<ValueType>::addGeometricalDependencies(
        int& testGeomDeps, int& trialGeomDeps) const
{
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
}

template <typename ValueType>
inline ValueType SingleLayerPotential3DKernel<ValueType>::evaluateAtPointPair(
        const arma::Col<ValueType>& testPoint,
        const arma::Col<ValueType>& trialPoint) const
{
    const int coordCount = testPoint.n_rows;

    ValueType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
    {
        ValueType diff = testPoint(coordIndex) - trialPoint(coordIndex);
        sum += diff * diff;
    }
    return 1. / (4. * M_PI * sqrt(sum));
}

template <typename ValueType>
void SingleLayerPotential3DKernel<ValueType>::evaluateAtPointPairs(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        arma::Cube<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;

#ifndef NDEBUG
    if (testPoints.n_rows != worldDimension() ||
            trialPoints.n_rows != worldDimension())
        throw std::invalid_argument("SingleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "3D coordinates required");
    if (testPoints.n_cols != trialPoints.n_cols)
        throw std::invalid_argument("SingleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "number of test and trial points must be equal");
#endif

    const int pointCount = testPoints.n_cols;
    result.set_size(1, 1, pointCount);
    for (int i = 0; i < pointCount; ++i)
        result(0, 0, i) = evaluateAtPointPair(testPoints.unsafe_col(i),
                                              trialPoints.unsafe_col(i));
}

template <typename ValueType>
void SingleLayerPotential3DKernel<ValueType>::evaluateOnGrid(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        Array4D<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;

#ifndef NDEBUG
    if (testPoints.n_rows != worldDimension() ||
            trialPoints.n_rows != worldDimension())
        throw std::invalid_argument("SingleLayerPotential3DKernel::evaluate(): "
                                    "3D coordinates required");
#endif

    const int testPointCount = testPoints.n_cols;
    const int trialPointCount = trialPoints.n_cols;
    result.set_size(1, testPointCount, 1, trialPointCount);
    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
            result(0, testIndex, 0, trialIndex) =
                    evaluateAtPointPair(testPoints.unsafe_col(testIndex),
                                        trialPoints.unsafe_col(trialIndex));
}

template<typename ValueType>
std::pair<const char*,int> SingleLayerPotential3DKernel<ValueType>::evaluateClCode () const
{
    return std::make_pair (single_layer_potential_3D_kernel_cl,
			   single_layer_potential_3D_kernel_cl_len);
}


#ifdef COMPILE_FOR_FLOAT
template class SingleLayerPotential3DKernel<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SingleLayerPotential3DKernel<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SingleLayerPotential3DKernel<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SingleLayerPotential3DKernel<std::complex<double> >;
#endif

} // namespace Fiber
