#include "standard_kernel_collection.hpp"

template <typename ValueType_>
class KernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef ... CoordinateType;

    int kernelCount() const;
    int kernelRowCount(int kernelIndex) const;
    int kernelColCount(int kernelIndex) const;

    void addGeometricalDependencies(int& testGeomDeps, int& trialGeomDeps) const;

    void evaluate(
            const GeometricDataSlice<CoordinateType>& testGeomData,
            const GeometricDataSlice<CoordinateType>& trialGeomData,
            Array3dCollectionSlice<ValueType>& result) const;
};

namespace Fiber
{

template <typename Functor>
void StandardKernelCollection<Functor>::evaluateAtPointPairs(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        Array3dCollection<ValueType>& result) const
{
#ifndef NDEBUG
    checkDimensions(testGeomData);
    checkDimensions(trialGeomData);
#endif
    assert(testGeomData.pointCount() == trialGeomData.pointCount());

    const int pointCount = testGeomData.pointCount();
    for (int i = 0; i < m_functor.kernelCount(); ++i)
        result.resizeElement(i,
                             m_functor.kernelRowCount(i),
                             m_functor.kernelColCount(i),
                             pointCount);

    for (int i = 0; i < pointCount; ++i)
        m_functor.evaluate(testGeomData.slice(i), trialGeomData.slice(i),
                           result.slice(i));
}

template <typename ValueType>
void Laplace3dSingleLayerPotentialKernel<ValueType>::evaluateOnGrid(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        Array4dCollection<ValueType>& result) const
{
#ifndef NDEBUG
    checkDimensions(testGeomData);
    checkDimensions(trialGeomData);
#endif
    assert(testGeomData.pointCount() == trialGeomData.pointCount());

    const int testPointCount = testGeomData.pointCount();
    const int trialPointCount = testGeomData.pointCount();
    // NOTE: this ordering doesn't match the one that is currently used
    for (int i = 0; i < m_functor.kernelCount(); ++i)
        result.resizeElement(i,
                             m_functor.kernelRowCount(i),
                             m_functor.kernelColCount(i),
                             testPointCount,
                             trialPointCount);

    for (int trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testPointCount; ++testIndex)
        m_functor.evaluate(testGeomData.slice(testIndex),
                           trialGeomData.slice(trialIndex),
                           result.slice(testIndex, trialIndex));
}

} // namespace Fiber
