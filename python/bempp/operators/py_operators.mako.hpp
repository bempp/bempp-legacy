
#ifndef BEMPP_OPERATORS_HPP
#define BEMPP_OPERATORS_HPP

#include "bempp/assembly/py_boundary_operator_variants.hpp"
#include "bempp/common/types.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/assembly/identity_operator.hpp"
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
    typedef shared_ptr<Space<BasisFunctionType>> space_t;

    space_t my_domain = boost::get<space_t>(domain.variants());
    space_t my_range = boost::get<space_t>(range.variants());
    space_t my_dual_to_range = boost::get<space_t>(dual_to_range.variants());

    return identityOperator<BasisFunctionType,ResultType>(
            parameterList,
            my_domain,my_range,my_dual_to_range,
            label,symmetry);
}

        
}
#endif
