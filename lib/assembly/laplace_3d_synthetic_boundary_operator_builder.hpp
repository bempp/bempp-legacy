#ifndef bempp_laplace_3d_synthetic_boundary_operator_builder_hpp
#define bempp_laplace_3d_synthetic_boundary_operator_builder_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <string>

namespace Bempp {

template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSyntheticBoundaryOperator(
    BoundaryOperator<BasisFunctionType, ResultType> (*constructor)(
        const shared_ptr<const Context<BasisFunctionType, ResultType>>
            & /*context*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*domain*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*range*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*dualToRange*/,
        const std::string & /*label*/, int /*symmetry*/),
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    std::string label, int internalSymmetry, int maximumSyntheseSymmetry);

} // namespace Bempp

#endif
