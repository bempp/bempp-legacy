#ifndef bempp_modified_helmholtz_3d_synthetic_boundary_operator_builder_hpp
#define bempp_modified_helmholtz_3d_synthetic_boundary_operator_builder_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <string>

namespace Bempp {

template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dSyntheticBoundaryOperator(
    BoundaryOperator<BasisFunctionType, ResultType>(*constructor)(
        const shared_ptr<
            const Context<BasisFunctionType, ResultType>> & /*context*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*domain*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*range*/,
        const shared_ptr<const Space<BasisFunctionType>> & /*dualToRange*/,
        KernelType /*waveNumber*/, const std::string & /*label*/,
        int /*symmetry*/, bool /*useInterpolation*/,
        int /*interpPtsPerWavelength*/),
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber, std::string label, int internalSymmetry,
    bool useInterpolation, int interpPtsPerWavelength,
    int maximumSyntheseSymmetry);

} // namespace Bempp

#endif
