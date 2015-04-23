<%
op_types = ['3dSingleLayerBoundaryOperator',
            '3dDoubleLayerBoundaryOperator',
            '3dAdjointDoubleLayerBoundaryOperator',
            '3dHypersingularBoundaryOperator']

%>

#ifndef BEMPP_OPERATORS_HPP
#define BEMPP_OPERATORS_HPP

#include "bempp/assembly/py_boundary_operator_variants.hpp"
#include "bempp/common/eigen_support.hpp"
#include "bempp/common/types.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_beltrami_3d_operator.hpp"
#include "bempp/assembly/maxwell_3d_identity_operator.hpp"
#include "bempp/assembly/maxwell_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/maxwell_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/maxwell_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_double_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "bempp/assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "bempp/assembly/discrete_boundary_operator.hpp"
#include "bempp/assembly/assembled_potential_operator.hpp"

#include "bempp/assembly/modified_helmholtz_3d_double_layer_potential_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_single_layer_potential_operator.hpp"

#include "bempp/assembly/helmholtz_3d_far_field_double_layer_potential_operator.hpp"
#include "bempp/assembly/helmholtz_3d_far_field_single_layer_potential_operator.hpp"

#include "bempp/assembly/maxwell_3d_far_field_double_layer_potential_operator.hpp"
#include "bempp/assembly/maxwell_3d_far_field_single_layer_potential_operator.hpp"


#include "boost/variant/get.hpp"
#include "bempp/common/shared_ptr.hpp"
#include "bempp/space/space.hpp"
#include "bempp/fiber/scalar_traits.hpp"



namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_identityOperator(
        const ParameterList& parameterList,
        const SpaceVariants& domain,
        const SpaceVariants& range,
        const SpaceVariants& dual_to_range,
        const std::string& label,
        int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);

    return identityOperator<BasisFunctionType,ResultType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,
            label,symmetry);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_laplaceBeltramiOperator(
        const ParameterList& parameterList,
        const SpaceVariants& domain,
        const SpaceVariants& range,
        const SpaceVariants& dual_to_range,
        const std::string& label,
        int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);

    return laplaceBeltrami3dOperator<BasisFunctionType,ResultType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,
            label,symmetry);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_laplace3dDoubleLayerBoundaryOperator(
                                                          const ParameterList& parameterList,
                                                          const SpaceVariants& domain,
                                                          const SpaceVariants& range,
                                                          const SpaceVariants& dual_to_range,
                                                          const std::string& label,
                                                          int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;
    
    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);
    
    return laplace3dDoubleLayerBoundaryOperator<BasisFunctionType,ResultType>(
                                                                              parameterList,
                                                                              my_domain,my_range,my_dual_to_range,
                                                                              label,symmetry);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_laplace3dSingleLayerBoundaryOperator(
                                                          const ParameterList& parameterList,
                                                          const SpaceVariants& domain,
                                                          const SpaceVariants& range,
                                                          const SpaceVariants& dual_to_range,
                                                          const std::string& label,
                                                          int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;
    
    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);
    
    return laplace3dSingleLayerBoundaryOperator<BasisFunctionType,ResultType>(
                                                                              parameterList,
                                                                              my_domain,my_range,my_dual_to_range,
                                                                              label,symmetry);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_laplace3dHypersingularBoundaryOperator(
                                                          const ParameterList& parameterList,
                                                          const SpaceVariants& domain,
                                                          const SpaceVariants& range,
                                                          const SpaceVariants& dual_to_range,
                                                          const std::string& label,
                                                          int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;
    
    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);
    
    return laplace3dHypersingularBoundaryOperator<BasisFunctionType,ResultType>(
                                                                              parameterList,
                                                                              my_domain,my_range,my_dual_to_range,
                                                                              label,symmetry);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_laplace3dAdjointDoubleLayerBoundaryOperator(
                                                          const ParameterList& parameterList,
                                                          const SpaceVariants& domain,
                                                          const SpaceVariants& range,
                                                          const SpaceVariants& dual_to_range,
                                                          const std::string& label,
                                                          int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;
    
    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);
    
    return laplace3dAdjointDoubleLayerBoundaryOperator<BasisFunctionType,ResultType>(
                                                                              parameterList,
                                                                              my_domain,my_range,my_dual_to_range,
                                                                              label,symmetry);
}

% for op in op_types:
template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_modifiedHelmholtz${op}(
                                       const ParameterList& parameterList,
                                       const SpaceVariants& domain,
                                       const SpaceVariants& range,
                                       const SpaceVariants& dual_to_range,
                                       ResultType waveNumber,
                                       const std::string& label,
                                       int symmetry)
{

    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);
    
    return modifiedHelmholtz${op}<BasisFunctionType,ResultType,ResultType>(
                                                                              parameterList,
                                                                              my_domain,my_range,my_dual_to_range,
                                                                              waveNumber,
                                                                              label,symmetry);

}
% endfor

static inline shared_ptr<const DiscreteBoundaryOperator<double>> 
py_laplace_single_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,double>> op( 
                new Laplace3dSingleLayerPotentialOperator<double,double>());
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}


static inline shared_ptr<const DiscreteBoundaryOperator<double>> 
py_laplace_double_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,double>> op( 
                new Laplace3dDoubleLayerPotentialOperator<double,double>());
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}

static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_modified_helmholtz_single_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new ModifiedHelmholtz3dSingleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}

static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_modified_helmholtz_double_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new ModifiedHelmholtz3dDoubleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOpVariants c_maxwellIdentityOperator(
        const ParameterList& parameterList,
        const SpaceVariants& domain,
        const SpaceVariants& range,
        const SpaceVariants& dual_to_range,
        const std::string& label,
        int symmetry)
{
    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);

    return maxwell3dIdentityOperator<BasisFunctionType,ResultType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,
            label,symmetry);
}

template <typename BasisFunctionType>
BoundaryOpVariants c_maxwellSingleLayerOperator(
        const ParameterList& parameterList,
        const SpaceVariants& domain,
        const SpaceVariants& range,
        const SpaceVariants& dual_to_range,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label,
        int symmetry)
{


    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);

    return maxwell3dSingleLayerBoundaryOperator<BasisFunctionType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,waveNumber,
            label,symmetry);

}

template <typename BasisFunctionType>
BoundaryOpVariants c_maxwellDoubleLayerOperator(
        const ParameterList& parameterList,
        const SpaceVariants& domain,
        const SpaceVariants& range,
        const SpaceVariants& dual_to_range,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label,
        int symmetry)
{


    typedef shared_ptr<const Space<BasisFunctionType>> space_t;

    space_t my_domain = _py_get_space_ptr<BasisFunctionType>(domain);
    space_t my_range = _py_get_space_ptr<BasisFunctionType>(range);
    space_t my_dual_to_range = _py_get_space_ptr<BasisFunctionType>(dual_to_range);

    return maxwell3dDoubleLayerBoundaryOperator<BasisFunctionType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,waveNumber,
            label,symmetry);

}


static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_maxwell_double_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Maxwell3dDoubleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}


static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_maxwell_single_layer_potential_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Maxwell3dSingleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}

static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_helmholtz_single_layer_far_field_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Helmholtz3dFarFieldSingleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}
static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_helmholtz_double_layer_far_field_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Helmholtz3dFarFieldDoubleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}
static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_maxwell_single_layer_far_field_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Maxwell3dFarFieldSingleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}
static inline shared_ptr<const DiscreteBoundaryOperator<std::complex<double>>> 
py_maxwell_double_layer_far_field_discrete_operator(
        const SpaceVariants& space, 
        const Matrix<double>& evaluationPoints,
        std::complex<double> waveNumber,
        const ParameterList& parameterList){

    typedef shared_ptr<const Space<double>> space_t;

    space_t my_space = _py_get_space_ptr<double>(space);
    shared_ptr<Matrix<double>> point_ptr(
            new Matrix<double>(evaluationPoints));

    shared_ptr<PotentialOperator<double,std::complex<double>>> op( 
                new Maxwell3dFarFieldDoubleLayerPotentialOperator<double>(
                    waveNumber));
    
    return op->assemble(my_space,
            point_ptr,
            parameterList).discreteOperator();
}
}
#endif
