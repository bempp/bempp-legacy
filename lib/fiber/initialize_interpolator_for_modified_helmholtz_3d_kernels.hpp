#ifndef bempp_initialize_interpolator_for_modified_helmholtz_3d_kernels_hpp
#define bempp_initialize_interpolator_for_modified_helmholtz_3d_kernels_hpp

#include "../common/common.hpp"
#include "hermite_interpolator.hpp"

namespace Fiber {

template <typename ValueType>
void initializeInterpolatorForModifiedHelmholtz3dKernels(
    ValueType waveNumber, typename ScalarTraits<ValueType>::RealType maxDist,
    int interpPtsPerWavelength, HermiteInterpolator<ValueType> &interpolator);

} // namespace Fiber

#endif
