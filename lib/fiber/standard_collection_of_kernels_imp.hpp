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
        size_t& testGeomDeps, size_t& trialGeomDeps) const
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

    const size_t pointCount = testGeomData.pointCount();
    const size_t kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           pointCount);

    for (size_t p = 0; p < pointCount; ++p)
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
    const size_t testPointCount = testGeomData.pointCount();
    const size_t trialPointCount = trialGeomData.pointCount();
    const size_t kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           testPointCount,
                           trialPointCount);

    for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex)
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
