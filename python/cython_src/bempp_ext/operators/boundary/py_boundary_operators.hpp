#ifndef bempp_py_boundary_operators_hpp
#define bempp_py_boundary_operators_hpp


#include "bempp/common/shared_ptr.hpp"
#include "bempp/assembly/abstract_boundary_operator.hpp"
#include "bempp/assembly/elementary_integral_operator.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_beltrami_3d_operator.hpp"
#include "bempp/assembly/maxwell_3d_identity_operator.hpp"
#include "bempp/assembly/maxwell_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/maxwell_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"

namespace Bempp {


inline shared_ptr<const ElementaryIntegralOperator<double,double,double>>
laplace_single_layer(const ParameterList& parameterList,
        const shared_ptr<const Space<double>>& domain,
        const shared_ptr<const Space<double>>& range,
        const shared_ptr<const Space<double>>& dualToRange,
        const std::string& label, int symmetry)
{
    return static_pointer_cast<const ElementaryIntegralOperator<double,double,double>>(
            laplace3dSingleLayerBoundaryOperator<double,double>(
                parameterList,domain,range,dualToRange,label,symmetry).
            abstractOperator());
}

inline shared_ptr<const ElementaryIntegralOperator<double,double,double>>
laplace_double_layer(const ParameterList& parameterList,
        const shared_ptr<const Space<double>>& domain,
        const shared_ptr<const Space<double>>& range,
        const shared_ptr<const Space<double>>& dualToRange,
        const std::string& label, int symmetry)
{
    return static_pointer_cast<const ElementaryIntegralOperator<double,double,double>>(
            laplace3dDoubleLayerBoundaryOperator<double,double>(
                parameterList,domain,range,dualToRange,label,symmetry).
            abstractOperator());
}

inline shared_ptr<const ElementaryIntegralOperator<double,double,double>>
laplace_adjoint_double_layer(const ParameterList& parameterList,
        const shared_ptr<const Space<double>>& domain,
        const shared_ptr<const Space<double>>& range,
        const shared_ptr<const Space<double>>& dualToRange,
        const std::string& label, int symmetry)
{
    return static_pointer_cast<const ElementaryIntegralOperator<double,double,double>>(
            laplace3dAdjointDoubleLayerBoundaryOperator<double,double>(
                parameterList,domain,range,dualToRange,label,symmetry).
            abstractOperator());
}

inline shared_ptr<const ElementaryIntegralOperator<double,double,double>>
laplace_hypersingular(const ParameterList& parameterList,
        const shared_ptr<const Space<double>>& domain,
        const shared_ptr<const Space<double>>& range,
        const shared_ptr<const Space<double>>& dualToRange,
        const std::string& label, int symmetry)
{
    return static_pointer_cast<const ElementaryIntegralOperator<double,double,double>>(
            laplace3dHypersingularBoundaryOperator<double,double>(
                parameterList,domain,range,dualToRange,label,symmetry).
            abstractOperator());
}
}


#endif
