#ifndef fiber_standard_collection_of_kernels_imp_hpp
#define fiber_standard_collection_of_kernels_imp_hpp

#include "standard_collection_of_kernels.hpp"

#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "geometrical_data.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename Functor>
void StandardCollectionOfKernels<Functor>::addGeometricalDependencies(
        int& testGeomDeps, int& trialGeomDeps) const
{
    m_functor.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
}

template <typename Functor>
void StandardCollectionOfKernels<Functor>::evaluateAtPointPairs(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<ValueType>& result) const
{
    assert(testGeomData.pointCount() == trialGeomData.pointCount());

    const int pointCount = testGeomData.pointCount();
    const int kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (int k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           pointCount);

    for (int p = 0; p < pointCount; ++p)
        m_functor.evaluate(testGeomData.const_slice(p),
                           trialGeomData.const_slice(p),
                           result.slice(p).self());
}

template <typename ValueType>
void StandardCollectionOfKernels<ValueType>::evaluateOnGrid(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf4dArrays<ValueType>& result) const
{
    const int testPointCount = testGeomData.pointCount();
    const int trialPointCount = trialGeomData.pointCount();
    const int kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (int k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           testPointCount,
                           trialPointCount);

    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
            m_functor.evaluate(testGeomData.const_slice(testIndex),
                               trialGeomData.const_slice(trialIndex),
                               result.slice(testIndex, trialIndex).self());
}

template <typename ValueType>
std::pair<const char*, int>
StandardCollectionOfKernels<ValueType>::evaluateClCode() const {
    throw std::runtime_error("StandardCollectionOfKernels::evaluateClCode(): "
                             "not implemented yet");
}

} // namespace Fiber

#endif
