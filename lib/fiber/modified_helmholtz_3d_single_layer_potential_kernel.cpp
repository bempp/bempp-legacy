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

#include "modified_helmholtz_3d_single_layer_potential_kernel.hpp"

#include "explicit_instantiation.hpp"

#include <armadillo>
#include <cmath>

#include "geometrical_data.hpp"
//#include "CL/modified_helmholtz_3d_single_layer_potential_kernel.cl.str"

namespace Fiber
{

template <typename ValueType>
ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::
ModifiedHelmholtz3dSingleLayerPotentialKernel(ValueType k)
{
    setWaveNumber(k);
}

template <typename ValueType>
void ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::addGeometricalDependencies(
        int& testGeomDeps, int& trialGeomDeps) const
{
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
}

template <typename ValueType>
inline ValueType ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::evaluateAtPointPair(
        const arma::Col<CoordinateType>& testPoint,
        const arma::Col<CoordinateType>& trialPoint) const
{
    const int coordCount = testPoint.n_rows;

    CoordinateType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
    {
        CoordinateType diff = testPoint(coordIndex) - trialPoint(coordIndex);
        sum += diff * diff;
    }
    CoordinateType distance = sqrt(sum);
    return static_cast<ValueType>(1.0 / (4.0*M_PI)) / distance * exp(-m_waveNumber*distance);
}

template <typename ValueType>
void ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::evaluateAtPointPairs(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        arma::Cube<ValueType>& result) const
{
    const arma::Mat<CoordinateType>& testPoints = testGeomData.globals;
    const arma::Mat<CoordinateType>& trialPoints = trialGeomData.globals;

#ifndef NDEBUG
    if (testPoints.n_rows != worldDimension() ||
            trialPoints.n_rows != worldDimension())
        throw std::invalid_argument("Laplace3dSingleLayerPotentialKernel::evaluateAtPointPairs(): "
                                    "3D coordinates required");
    if (testPoints.n_cols != trialPoints.n_cols)
        throw std::invalid_argument("Laplace3dSingleLayerPotentialKernel::evaluateAtPointPairs(): "
                                    "number of test and trial points must be equal");
#endif

    const int pointCount = testPoints.n_cols;
    result.set_size(1, 1, pointCount);
    for (int i = 0; i < pointCount; ++i)
        result(0, 0, i) = evaluateAtPointPair(testPoints.unsafe_col(i),
                                              trialPoints.unsafe_col(i));
}

template <typename ValueType>
void ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::evaluateOnGrid(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        Array4d<ValueType>& result) const
{
    const arma::Mat<CoordinateType>& testPoints = testGeomData.globals;
    const arma::Mat<CoordinateType>& trialPoints = trialGeomData.globals;

#ifndef NDEBUG
    if (testPoints.n_rows != worldDimension() ||
            trialPoints.n_rows != worldDimension())
        throw std::invalid_argument("Laplace3dSingleLayerPotentialKernel::evaluate(): "
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
std::pair<const char*,int> ModifiedHelmholtz3dSingleLayerPotentialKernel<ValueType>::evaluateClCode () const
{
    return std::make_pair ("", 0);  // TODO

  //    return std::make_pair (modified_helmholtz_3d_single_layer_potential_kernel_cl,
  //			   modified_helmholtz_3d_single_layer_potential_kernel_cl_len);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_KERNEL(ModifiedHelmholtz3dSingleLayerPotentialKernel);

} // namespace Fiber
