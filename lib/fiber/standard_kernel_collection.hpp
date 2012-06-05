#ifndef fiber_standard_kernel_collection_hpp
#define fiber_standard_kernel_collection_hpp

namespace Fiber
{

template <typename Functor>
class StandardKernelCollection :
        public KernelCollection<typename Functor::ValueType>
{
public:
    typedef typename Functor::ValueType ValueType;
    virtual void evaluateAtPointPairs(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            Array3dCollection<ValueType>& result) const;

    virtual void evaluateOnGrid(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            Array4dCollection<ValueType>& result) const;
};

} // namespace Fiber

#endif
