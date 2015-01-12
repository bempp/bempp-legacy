
#ifndef BEMPP_OPERATORS_HPP
#define BEMPP_OPERATORS_HPP

#include "bempp/assembly/py_boundary_operator_variants.hpp"
#include "bempp/common/types.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "boost/variant/get.hpp"
#include "bempp/common/shared_ptr.hpp"
#include "bempp/space/space.hpp"


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


}
#endif
