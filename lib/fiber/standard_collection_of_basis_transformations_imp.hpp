#include "standard_collection_of_basis_transformations.hpp"

namespace Fiber
{

template <typename Functor>
void StandardCollectionOfBasisTransformations<Functor>::addDependencies(
        int& basisDeps, int& geomDeps) const
{
    m_functor.addDependencies(basisDeps, geomDeps);
}

template <typename Functor>
int StandardCollectionOfBasisTransformations<Functor>::
transformationCount() const
{
    return m_functor.transformationCount();
}

template <typename Functor>
int StandardCollectionOfBasisTransformations<Functor>::
argumentDimension() const
{
    return m_functor.argumentDimension();
}

template <typename Functor>
int StandardCollectionOfBasisTransformations<Functor>::
resultDimension(int transformationIndex) const
{
    return m_functor.resultDimension(transformationIndex);
}

template <typename Functor>
template <typename ValueType>
void StandardCollectionOfBasisTransformations<Functor>::evaluateImpl(
        const BasisData<ValueType>& basisData,
        const GeometricalData<CoordinateType>& geomData,
        CollectionOf3dArrays<ValueType>& result) const
{
    const int pointCount = basisData.pointCount();
    const int functionCount = basisData.functionCount();
    const int transformationCount = m_functor.transformationCount();
    result.set_size(transformationCount);
    for (int t = 0; t < transformationCount; ++t)
        result[t].set_size(m_functor.resultDimension(t),
                           functionCount,
                           pointCount);

    for (int p = 0; p < pointCount; ++p)
        for (int f = 0; f < functionCount; ++f)
            m_functor.evaluate(basisData.const_slice(f, p),
                               geomData.const_slice(p),
                               result.slice(f, p).self());
}

template <typename Functor>
void StandardCollectionOfBasisTransformations<Functor>::evaluateImplReal(
                const BasisData<CoordinateType>& basisData,
                const GeometricalData<CoordinateType>& geomData,
                CollectionOf3dArrays<CoordinateType>& result) const
{
    evaluateImpl(basisData, geomData, result);
}

template <typename Functor>
void StandardCollectionOfBasisTransformations<Functor>::evaluateImplComplex(
                const BasisData<ComplexType>& basisData,
                const GeometricalData<CoordinateType>& geomData,
                CollectionOf3dArrays<ComplexType>& result) const
{
    evaluateImpl(basisData, geomData, result);
}

} // namespace Fiber
