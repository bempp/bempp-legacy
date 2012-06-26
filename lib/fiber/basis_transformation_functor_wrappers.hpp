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

#ifndef fiber_functor_wrappers_for_default_collection_of_basis_transformations_hpp
#define fiber_functor_wrappers_for_default_collection_of_basis_transformations_hpp

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"

namespace Fiber
{

// ElementaryFunctor must be default-constructible
template <typename ElementaryFunctor>
class ElementaryBasisTransformationFunctorWrapper
{
public:
    typedef typename ElementaryFunctor::CoordinateType CoordinateType;

    size_t transformationCount() const {
        return 1;
    }

    int argumentDimension() const {
        return m_functor.argumentDimension();
    }

    int resultDimension(size_t transformationIndex) const {
        assert(transformationIndex == 0);
        return m_functor.resultDimension();
    }

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        m_functor.addDependencies(basisDeps, geomDeps);
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            CollectionOf1dSlicesOf3dArrays<ValueType>& result) const {
        m_functor.evaluate(basisData, geomData, result[0].self());
    }

private:
    ElementaryFunctor m_functor;
};

// ElementaryFunctors must be default-constructible
template <typename ElementaryFunctor0, typename ElementaryFunctor1>
class ElementaryBasisTransformationFunctorPairWrapper
{
public:
    typedef typename ElementaryFunctor0::CoordinateType CoordinateType;
    // TODO: check that ElementaryFunctor1 has the same CoordinateType

    int transformationCount() const {
        return 2;
    }

    int argumentDimension() const {
        int result = m_functor0.argumentDimension();
        assert(result == m_functor1.argumentDimension());
        return result;
    }

    int resultDimension(int transformationIndex) const {
        switch (transformationIndex) {
        case 0: return m_functor0.resultDimension();
        case 1: return m_functor1.resultDimension();
        default: throw std::invalid_argument(
                        "BasisFunctionTransformationElementaryFunctorPairWrapper::"
                        "resultDimension(): invalid transformation index");
        }
    }

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        m_functor0.addDependencies(basisDeps, geomDeps);
        m_functor1.addDependencies(basisDeps, geomDeps);
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            CollectionOf1dSlicesOf3dArrays<ValueType>& result) const {
        m_functor0.evaluate(basisData, geomData, result[0].self());
        m_functor1.evaluate(basisData, geomData, result[1].self());
    }

private:
    ElementaryFunctor0 m_functor0;
    ElementaryFunctor1 m_functor1;
};

} // namespace Fiber

#endif
