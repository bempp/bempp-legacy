#ifndef fiber_standard_collection_of_basis_transformations_hpp
#define fiber_standard_collection_of_basis_transformations_hpp

#include "collection_of_basis_transformations.hpp"

namespace Fiber
{

/*
template <typename CoordinateType_>
class BasisFunctionTransformationFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int transformationCount() const;
    int transformationDim(int transformationIndex) const;

    void addDependencies(int& basisDeps, int& geomDeps) const;

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            Array3dNonconstSlice1dCollection<ValueType>& result) const;
};
*/

template <typename Functor>
class StandardCollectionOfBasisTransformations :
        public CollectionOfBasisTransformations<typename Functor::CoordinateType>
{
    typedef CollectionOfBasisTransformations<typename Functor::CoordinateType>
    Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::ComplexType ComplexType;

    StandardCollectionOfBasisTransformations(const Functor& functor) :
        m_functor(functor)
    {}

    virtual int transformationCount() const;
    virtual int argumentDimension() const;
    virtual int resultDimension(int transformationIndex) const;
    virtual void addDependencies(size_t& basisDeps, size_t& geomDeps) const;

private:
    template <typename ValueType>
    void evaluateImpl(
            const BasisData<ValueType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<ValueType>& result) const;

    virtual void evaluateImplReal(
            const BasisData<CoordinateType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<CoordinateType>& result) const;

    virtual void evaluateImplComplex(
            const BasisData<ComplexType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<ComplexType>& result) const;

private:
    Functor m_functor;
};

} // namespace Fiber

#include "standard_collection_of_basis_transformations_imp.hpp"

#endif
