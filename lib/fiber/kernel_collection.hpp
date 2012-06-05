#ifndef fiber_kernel_collection_hpp
#define fiber_kernel_collection_hpp

namespace Fiber
{

template <typename ValueType>
class KernelCollection
{
    virtual void evaluateAtPointPairs(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            Array3dCollection<ValueType>& result) const = 0;

    virtual void evaluateOnGrid(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            Array4dCollection<ValueType>& result) const = 0;

    /**
     * \brief Returns an OpenCL code snippet for kernel evaluation as a string.
     * \note The code snippet must provide device functions devKernevalGrid and
     *    devKernevalPair (for an implementation example, see
     *    CL/laplace_3d_single_layer_potential_kernel.cl)
     * \note The required data must have been pushed to device memory before
     *    invocation.
     */
    virtual std::pair<const char*,int> evaluateClCode() const = 0;
};

} // namespace Fiber

#endif
