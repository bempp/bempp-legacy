#include "single_layer_potential_3d_kernel.hpp"

#include <armadillo>
#include <cmath>

#include "geometrical_data.hpp"

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
    const int coordCount = testPoints.n_rows;
    result.set_size(1, 1, testPointCount,  trialPointCount);
    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
            result(0, 0, testIndex, trialIndex) =
                    evaluateAtPointPair(testPoints.unsafe_col(testIndex),
                                        trialPoints.unsafe_col(trialIndex));
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
