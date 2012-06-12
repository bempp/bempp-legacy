#ifndef fiber_standard_collection_of_kernels_hpp
#define fiber_standard_collection_of_kernels_hpp

#include "collection_of_kernels.hpp"

namespace Fiber
{

/*
template <typename ValueType_>
class KernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int kernelCount() const;
    int kernelRowCount(int kernelIndex) const;
    int kernelColCount(int kernelIndex) const;

    void addGeometricalDependencies(int& testGeomDeps, int& trialGeomDeps) const;

    template <template <typename T> class ArrayNdNonconstSlice2dCollection>
    void evaluate(
            const GeometricDataSlice<CoordinateType>& testGeomData,
            const GeometricDataSlice<CoordinateType>& trialGeomData,
            ArrayNdNonconstSlice2dCollection<ValueType>& result) const;
};
*/

template <typename Functor>
class StandardCollectionOfKernels :
        public CollectionOfKernels<typename Functor::ValueType>
{
    typedef CollectionOfKernels<typename Functor::ValueType> Base;
public:
    typedef typename Base::ValueType ValueType;
    typedef typename Base::CoordinateType CoordinateType;

    StandardCollectionOfKernels(const Functor& functor) :
        m_functor(functor)
    {}

    const Functor& functor() const {
        return m_functor;
    }

    Functor& functor() {
        return m_functor;
    }

    virtual void addGeometricalDependencies(
            size_t& testGeomDeps, size_t& trialGeomDeps) const;

    virtual void evaluateAtPointPairs(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf3dArrays<ValueType>& result) const;

    virtual void evaluateOnGrid(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf4dArrays<ValueType>& result) const;

    virtual std::pair<const char*, int> evaluateClCode() const;

private:
    Functor m_functor;
};

} // namespace Fiber

#include "standard_collection_of_kernels_imp.hpp"

#endif
