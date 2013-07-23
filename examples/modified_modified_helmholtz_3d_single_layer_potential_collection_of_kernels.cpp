#include "modified_modified_helmholtz_3d_single_layer_potential_collection_of_kernels.hpp"

#include "common/complex_aux.hpp"
#include "fiber/collection_of_3d_arrays.hpp"
#include "fiber/collection_of_4d_arrays.hpp"
#include "fiber/explicit_instantiation.hpp"
#include "modified_geometrical_data.hpp"
#include "modified_kernel_values.hpp"

#include "common/armadillo_fwd.hpp"
#include <stdexcept>

template <typename ValueType>
inline ValueType myexpm(const ValueType& x)
{
    return std::exp(-x);
}

template <typename ValueType>
inline std::complex<ValueType> myexpm(const std::complex<ValueType>& x)
{
    return std::exp(-x.real()) * 
        std::complex<ValueType>(cos(x.imag()), -sin(x.imag()));
}

namespace Fiber
{

template <typename ValueType>
void ModifiedModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels<ValueType>::evaluateOnGrid(
    const ModifiedGeometricalData<CoordinateType>& /*__restrict*/ testGeomData,
        const ModifiedGeometricalData<CoordinateType>& /*__restrict*/ trialGeomData,
        ModifiedKernelValues<ValueType>& /*__restrict*/ result) const
{
    const size_t testPointCount = testGeomData.pointCount;
    const size_t trialPointCount = trialGeomData.pointCount;

    // ValueType* __restrict **  values_ = result.values;
    ValueType* __restrict  values00_ = result.values[0][0];
    const CoordinateType* __restrict testGlobals[] = {testGeomData.globals[0], testGeomData.globals[1], testGeomData.globals[2]};
    const CoordinateType* __restrict trialGlobals[] = {trialGeomData.globals[0], trialGeomData.globals[1], trialGeomData.globals[2]};

    // const CoordinateType** __restrict trialGlobals = trialGeomData.globals;
    // std::cout << (size_t)result_ % 16 << std::endl;
    // std::cout << (size_t)testGlobals % 16 << std::endl;
    // std::cout << (size_t)trialGlobals % 16 << std::endl;
    // __assume_aligned(result.values, 16);
    // __assume_aligned(testGeomData.globals, 16);
    // __assume_aligned(trialGeomData.globals, 16);

    const int coordCount = 3;
// #pragma vector aligned
//#pragma ivdep
    for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
//#pragma ivdep
       for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex) {
           CoordinateType sum = 0;
           for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
           {
               CoordinateType diff =
                   // testGeomData.globals[coordIndex][testIndex] -
                   // trialGeomData.globals[coordIndex][trialIndex];
                   testGlobals[coordIndex][testIndex] -
                   trialGlobals[coordIndex][trialIndex];
               sum += diff * diff;
           }
           CoordinateType distance = sqrt(sum);
           // result.values[0][0][testIndex + testPointCount * trialIndex] =
           // values_[0][0][testIndex + testPointCount * trialIndex] =
           values00_[testIndex + testPointCount * trialIndex] =
               static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance *
               myexpm(m_waveNumber * distance);
       }


}

template class ModifiedModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels<std::complex<double> >;

// FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(ModifiedModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels);

} // namespace Fiber
