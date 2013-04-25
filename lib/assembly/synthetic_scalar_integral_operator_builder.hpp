#ifndef bempp_synthetic_scalar_integral_operator_builder_hpp
#define bempp_synthetic_scalar_integral_operator_builder_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "symmetry.hpp"

#include <string>

namespace Bempp
{

template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator;

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
makeSyntheticScalarIntegralOperator(
        const BoundaryOperator<BasisFunctionType, ResultType>& internalOp,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const shared_ptr<const Space<BasisFunctionType> >& internalTrialSpace,
        const shared_ptr<const Space<BasisFunctionType> >& internalTestSpace,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

} // namespace Bempp

#endif
