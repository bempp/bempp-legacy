#ifndef bempp_vc_modified_helmholtz_3d_single_layer_potential_collection_of_kernels_hpp
#define bempp_vc_modified_helmholtz_3d_single_layer_potential_collection_of_kernels_hpp

#include "fiber/collection_of_4d_arrays.hpp"
#include "fiber/scalar_traits.hpp"
#include <Vc/Vc>

namespace Fiber
{

template <typename T> class VcGeometricalData;
template <typename T> class VcKernelValues;

template <typename ValueType>
class VcModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels
{
public:
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    explicit VcModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels(ValueType k, ValueType kImag=0.) : m_waveNumber(k), m_waveNumberImag(kImag) {}

void evaluateOnGrid(
    const VcGeometricalData<CoordinateType>& testGeomData,
    const VcGeometricalData<CoordinateType>& trialGeomData,
    CollectionOf4dArrays<Vc::Vector<ValueType> >& result,
    CollectionOf4dArrays<Vc::Vector<ValueType> >& resultImag) const;

private:
    CoordinateType m_waveNumber, m_waveNumberImag;
};

}

#endif
