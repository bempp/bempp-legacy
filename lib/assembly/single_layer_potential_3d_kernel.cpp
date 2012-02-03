#include "single_layer_potential_3d_kernel.hpp"

#include <armadillo>
#include <cmath>

namespace Bempp
{

template <typename ValueType>
void SingleLayerPotential3DKernel<ValueType>::evaluate(
        const arma::Mat<ctype>& testPoints,
        const arma::Mat<ctype>& trialPoints,
        const arma::Mat<ctype>& /* testNormals */,
        const arma::Mat<ctype>& /* trialNormals */,
        arma::Cube<ValueType>& result) const
{
#ifndef NDEBUG
    if (testPoints.n_rows != worldDimension() ||
            trialPoints.n_rows != worldDimension())
        throw std::invalid_argument("SingleLayerPotential3DKernel::evaluate(): "
                                    "3D coordinates required");
    if (testPoints.n_cols != trialPoints.n_cols)
        throw std::invalid_argument("SingleLayerPotential3DKernel::evaluate(): "
                                    "number of test and trial points must be equal");
#endif
    const int pointCount = testPoints.n_cols;
    const int coordCount = testPoints.n_rows;
    result.set_size(1, 1, pointCount);
    for (int i = 0; i < pointCount; ++i)
    {
        ValueType sum = 0;
        for (int j = 0; j < coordCount; ++j)
            sum += testPoints(j, i) * testPoints(j, i);
        result(1, 1, i) = sqrt(sum);
    }
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

} // namespace Bempp
