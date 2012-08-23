// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef fiber_default_collection_of_basis_transformations_hpp
#define fiber_default_collection_of_basis_transformations_hpp

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
class DefaultCollectionOfBasisTransformations :
        public CollectionOfBasisTransformations<typename Functor::CoordinateType>
{
    typedef CollectionOfBasisTransformations<typename Functor::CoordinateType>
    Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::ComplexType ComplexType;

    explicit DefaultCollectionOfBasisTransformations(const Functor& functor) :
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

#include "default_collection_of_basis_transformations_imp.hpp"

#endif
