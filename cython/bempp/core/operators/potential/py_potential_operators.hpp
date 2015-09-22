#ifndef bempp_core_py_potential_operators_hpp
#define bempp_core_py_potential_operators_hpp

#include "bempp/common/shared_ptr.hpp"
#include "bempp/common/eigen_support.hpp"
#include "bempp/assembly/assembled_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_potential_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_double_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_double_layer_potential_operator.hpp"

#include <complex>


namespace Bempp {


inline shared_ptr<const DiscreteBoundaryOperator<double>> laplace_single_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, double>> op(
                new Laplace3dSingleLayerPotentialOperator<double>());
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<double>> laplace_double_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, double>> op(
                new Laplace3dDoubleLayerPotentialOperator<double>());
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}


inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> modified_helmholtz_single_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new ModifiedHelmholtz3dSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> modified_helmholtz_double_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new ModifiedHelmholtz3dDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}


inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> maxwell_single_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> maxwell_double_layer_potential_discrete_operator(
        const shared_ptr<const Space<double>>& space,
        const Matrix<double>& evaluationPoints,
        const std::complex<double> waveNumber, 
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

}

#endif
