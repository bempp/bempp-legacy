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

#ifndef fiber_expression_hpp
#define fiber_expression_hpp

#include "scalar_traits.hpp"

#include <armadillo>

namespace Fiber {

template <typename ValueType> class BasisData;
template <typename CoordinateType> class GeometricalData;

template <typename CoordinateType>
class Expression
{
public:
    typedef typename ScalarTraits<CoordinateType>::ComplexType ComplexType;

    virtual ~Expression() {}

    virtual size_t domainDimension() const = 0;
    virtual size_t codomainDimension() const = 0;

    virtual void addDependencies(size_t& basisDeps, size_t& geomDeps) const = 0;

    /** \brief Evaluate expression at specified points.
     *
     *  \param[out] result
     *    Result. Array dimensions: (codomainDim, functionCount, pointCount). */
    void evaluate(const BasisData<CoordinateType>& basisData,
                  const GeometricalData<CoordinateType>& geomData,
                  arma::Cube<CoordinateType>& result) const {
        evaluateImplReal(basisData, geomData, result);
    }

    void evaluate(const BasisData<ComplexType>& basisData,
                  const GeometricalData<CoordinateType>& geomData,
                  arma::Cube<ComplexType>& result) const {
        evaluateImplComplex(basisData, geomData, result);
    }

private:
    virtual void evaluateImplReal(
            const BasisData<CoordinateType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            arma::Cube<CoordinateType>& result) const = 0;

    virtual void evaluateImplComplex(
            const BasisData<ComplexType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            arma::Cube<ComplexType>& result) const = 0;
};

} // namespace Fiber

#endif
