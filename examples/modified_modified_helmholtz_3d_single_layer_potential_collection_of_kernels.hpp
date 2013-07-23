#include "fiber/scalar_traits.hpp"

namespace Fiber
{

template <typename T> class ModifiedGeometricalData;
template <typename T> class ModifiedKernelValues;

template <typename ValueType>
class ModifiedModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels
{
public:
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    explicit ModifiedModifiedHelmholtz3dSingleLayerPotentialCollectionOfKernels(ValueType k) : m_waveNumber(k) {}

void evaluateOnGrid(
    const ModifiedGeometricalData<CoordinateType>& testGeomData,
    const ModifiedGeometricalData<CoordinateType>& trialGeomData,
    ModifiedKernelValues<ValueType>& result) const;

private:
    ValueType m_waveNumber;
};

}
