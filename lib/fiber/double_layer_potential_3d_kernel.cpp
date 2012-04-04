#include "double_layer_potential_3d_kernel.hpp"

#include <armadillo>
#include <cassert>
#include <cmath>

#include "geometrical_data.hpp"
#include "CL/double_layer_potential_3D_kernel.cl.str"

namespace Fiber
{

// Double potential: derivative wrt. trial normal

template <typename ValueType>
void DoubleLayerPotential3DKernel<ValueType>::addGeometricalDependencies(
        int& testGeomDeps, int& trialGeomDeps) const
{
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS | NORMALS;
}

template <typename ValueType>
inline ValueType DoubleLayerPotential3DKernel<ValueType>::evaluateAtPointPair(
        const arma::Col<ValueType>& testPoint,
        const arma::Col<ValueType>& trialPoint,
        const arma::Col<ValueType>& trialNormal) const
{
    const int coordCount = testPoint.n_rows;

    ValueType numeratorSum = 0., denominatorSum = 0.;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
    {
        ValueType diff = testPoint(coordIndex) - trialPoint(coordIndex);
        denominatorSum += diff * diff;
        numeratorSum += diff * trialNormal(coordIndex);
    }
    ValueType distance = sqrt(denominatorSum);
    return -numeratorSum / (4. * M_PI * distance * distance * distance);
}

template <typename ValueType>
void DoubleLayerPotential3DKernel<ValueType>::evaluateAtPointPairs(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        arma::Cube<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;
    const arma::Mat<ValueType>& trialNormals = trialGeomData.normals;

#ifndef NDEBUG
    const int worldDim = worldDimension();
    if (testPoints.n_rows != worldDim || trialPoints.n_rows != worldDim)
        throw std::invalid_argument("DoubleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "3D coordinates required");
    if (testPoints.n_cols != trialPoints.n_cols)
        throw std::invalid_argument("DoubleLayerPotential3DKernel::evaluateAtPointPairs(): "
                                    "number of test and trial points must be equal");
    assert(trialNormals.n_rows == worldDim);
    assert(trialNormals.n_cols == trialPoints.n_cols);
#endif

    const int pointCount = testPoints.n_cols;
    result.set_size(1, 1, pointCount);
    for (int i = 0; i < pointCount; ++i)
        result(0, 0, i) = evaluateAtPointPair(
                    testPoints.unsafe_col(i), trialPoints.unsafe_col(i),
                    trialNormals.unsafe_col(i));
}

template <typename ValueType>
void DoubleLayerPotential3DKernel<ValueType>::evaluateOnGrid(
        const GeometricalData<ValueType>& testGeomData,
        const GeometricalData<ValueType>& trialGeomData,
        Array4D<ValueType>& result) const
{
    const arma::Mat<ValueType>& testPoints = testGeomData.globals;
    const arma::Mat<ValueType>& trialPoints = trialGeomData.globals;
    const arma::Mat<ValueType>& trialNormals = trialGeomData.normals;

#ifndef NDEBUG
    const int worldDim = worldDimension();
    if (testPoints.n_rows != worldDim || trialPoints.n_rows != worldDim)
        throw std::invalid_argument("DoubleLayerPotential3DKernel::evaluate(): "
                                    "3D coordinates required");
    assert(trialNormals.n_rows == worldDim);
    assert(trialNormals.n_cols == trialPoints.n_cols);
#endif

    const int testPointCount = testPoints.n_cols;
    const int trialPointCount = trialPoints.n_cols;
    result.set_size(1, 1, testPointCount,  trialPointCount);
    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
            result(0, 0, testIndex, trialIndex) = evaluateAtPointPair(
                        testPoints.unsafe_col(testIndex),
                        trialPoints.unsafe_col(trialIndex),
                        trialNormals.unsafe_col(trialIndex));
}

template<typename ValueType>
std::pair<const char*,int> DoubleLayerPotential3DKernel<ValueType>::evaluateClCode () const
{
    return std::make_pair(double_layer_potential_3D_kernel_cl,
			  double_layer_potential_3D_kernel_cl_len);
}

#ifdef COMPILE_FOR_FLOAT
template class DoubleLayerPotential3DKernel<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class DoubleLayerPotential3DKernel<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class DoubleLayerPotential3DKernel<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class DoubleLayerPotential3DKernel<std::complex<double> >;
#endif

} // namespace Fiber
