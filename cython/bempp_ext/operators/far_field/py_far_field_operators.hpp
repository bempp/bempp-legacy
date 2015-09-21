#ifndef bempp_ext_py_far_field_operators_hpp
#define bempp_ext_py_far_field_operators_hpp

#include "bempp/common/shared_ptr.hpp"
#include "bempp/common/eigen_support.hpp"
#include "bempp/assembly/assembled_potential_operator.hpp"
#include "bempp/assembly/helmholtz_3d_far_field_single_layer_potential_operator.hpp"
#include "bempp/assembly/helmholtz_3d_far_field_double_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_far_field_single_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_far_field_double_layer_potential_operator.hpp"

#include <complex>


namespace Bempp {


inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> helmholtz_single_layer_far_field_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Helmholtz3dFarFieldSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> helmholtz_double_layer_far_field_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Helmholtz3dFarFieldDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> maxwell_single_layer_far_field_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dFarFieldSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> maxwell_double_layer_far_field_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dFarFieldDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

}

#endif
