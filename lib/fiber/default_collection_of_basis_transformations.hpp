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

/** \ingroup weak_form_elements
 *  \brief Default implementation of a collection of basis function transformations.

    This class implements the interface defined by
    CollectionOfBasisTransformations using a functor object to evaluate the
    transformations at specific points.

    \tparam Functor
      Type of the functor that will be passed to the constructor and used to
      evaluate a number of basis function transformations at individual points.

    The Functor class should provide the following interface:

    \code
    class TransformationCollectionFunctor
    {
    public:
        typedef ... CoordinateType; // should be either float or double

        // Return the number of transformations in the collection
        int transformationCount() const;
        // Return the dimension of vectors being the values of basis functions
        // to be transformed
        int argumentDimension() const;
        // Return the dimension of vectors produced by i'th transformation
        int resultDimension(int i) const;

        // Specify types of data required by the transformations (see the
        // documentation of CollectionOfBasisTransformations::addDependencies()
        // for details)
        void addDependencies(size_t& basisDeps, size_t& geomDeps) const;

        // Evaluate the transformation of a basis function whose values and/or
        // derivatives are provided in the basisData argument, at the point
        // whose geometrical data are provided in the geomData argument. The
        // j'th component of the vector being the value of i'th transformation
        // of the basis function should be written to result[i](j).
        // This function should accept ValueType equal to either
        // CoordinateType or std::complex<CoordinateType>
        template <typename ValueType>
        void evaluate(
                const ConstBasisDataSlice<ValueType>& basisData,
                const ConstGeometricalDataSlice<CoordinateType>& geomData,
                CollectionOf1dSlicesOf3dArrays<ValueType>& result) const;
    };
   \endcode
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
