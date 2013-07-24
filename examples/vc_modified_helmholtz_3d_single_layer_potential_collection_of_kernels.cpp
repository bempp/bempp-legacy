#include "vc_modified_helmholtz_3d_single_layer_potential_collection_of_kernels.hpp"

#include "common/complex_aux.hpp"
#include "fiber/explicit_instantiation.hpp"
#include "vc_geometrical_data.hpp"
#include "vc_multicomponent_kernel_values.hpp"

#include <Vc/Vc>
#include <stdexcept>

namespace Fiber
{

template <typename CoordinateType>
inline void evaluateHlmKernel(const Vc::Vector<CoordinateType>** testGlobalCoords,
                              const Vc::Vector<CoordinateType>** trialGlobalCoords,
                              CollectionOf2dSectionsOf4dArrays<Vc::Vector<CoordinateType> >& result,
                              CollectionOf2dSectionsOf4dArrays<Vc::Vector<CoordinateType> >& resultImag,
                              CoordinateType waveNumber,
                              CoordinateType waveNumberImag)
{
    const int coordCount = 3;

    Vc::Vector<CoordinateType> diff;
    Vc::Vector<CoordinateType> distance;
    Vc::Vector<CoordinateType> exponential, sine, cosine;
    Vc::Vector<CoordinateType> sum(Vc::Zero);    
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
        diff = *testGlobalCoords[coordIndex] - *trialGlobalCoords[coordIndex];
        sum = sum + diff * diff;
    }
    distance = Vc::sqrt(sum);
    exponential = Vc::exp(-waveNumberImag * distance);
    exponential *= (static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance);
    result[0](0, 0) = exponential;
    resultImag[0](0, 0) = exponential;            
    distance *= waveNumber; 
    Vc::sincos(distance, &sine, &cosine);
    result[0](0, 0) *= cosine;
    resultImag[0](0, 0) *= sine;
}

template <typename ValueType>
void VcModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels<ValueType>::evaluateOnGrid(
        const VcGeometricalData<CoordinateType>& testGeomData,
        const VcGeometricalData<CoordinateType>& trialGeomData,
        CollectionOf4dArrays<Vc::Vector<ValueType> >& result,
        CollectionOf4dArrays<Vc::Vector<ValueType> >& resultImag) const
{
    const int coordCount = 3;
    const size_t vsize = Vc::Vector<CoordinateType>::Size;
    const size_t testPointChunkCount = testGeomData.paddedPointCount / vsize;
    const size_t trialPointCount = trialGeomData.pointCount;
    const size_t paddedTrialPointCount = trialGeomData.paddedPointCount;

    Vc::Vector<CoordinateType> trialGlobalCoordsMem[coordCount];
    const Vc::Vector<CoordinateType>* testGlobalCoords[coordCount];
    const Vc::Vector<CoordinateType>* trialGlobalCoords[coordCount];
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
        trialGlobalCoords[coordIndex] = trialGlobalCoordsMem + coordIndex;
    
    for (size_t trialIndex = 0, resultIndex = 0; trialIndex < trialPointCount; ++trialIndex) {
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
            trialGlobalCoordsMem[coordIndex] = 
                trialGeomData.globals[trialIndex + coordIndex * paddedTrialPointCount];
        
        for (size_t testVIndex = 0; testVIndex < testPointChunkCount;
             ++testVIndex, ++resultIndex) {
            for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
                testGlobalCoords[coordIndex] = 
                    reinterpret_cast<const Vc::Vector<CoordinateType>*>(
                        testGeomData.globals.entries()) + 
                    testVIndex + coordIndex * testPointChunkCount;
            evaluateHlmKernel(testGlobalCoords, trialGlobalCoords,
                              result.section(testVIndex, trialIndex).self(),
                              resultImag.section(testVIndex, trialIndex).self(),
                              m_waveNumber, m_waveNumberImag);
        }
    }
}

template class VcModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels<double>;

} // namespace Fiber
