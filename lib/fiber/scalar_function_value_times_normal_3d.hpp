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

#ifndef fiber_scalar_function_value_times_normal_3d_hpp
#define fiber_scalar_function_value_times_normal_3d_hpp

#include "../common/common.hpp"

#include "expression.hpp"
#include "scalar_space_mapping.hpp"

#include "../common/armadillo_fwd.hpp"

namespace Fiber
{

/** \brief Expression defined as the vector normal to (3D) surface multiplied
 *  by a scalar function. */
template <typename CoordinateType>
class ScalarFunctionValueTimesNormal3d : public Expression<CoordinateType>
{
public:
    typedef typename Expression<CoordinateType>::ComplexType ComplexType;

    virtual size_t domainDimension() const {
        return 1;
    }

    virtual size_t codomainDimension() const {
        return 3;
    }

    virtual void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        ScalarSpaceMapping<CoordinateType>::
                addShapeFunctionDependencies(basisDeps, geomDeps);
        geomDeps |= NORMALS;
    }

private:
    template <typename ValueType>
    void evaluateImpl(const BasisData<ValueType>& basisData,
                      const GeometricalData<CoordinateType>& geomData,
                      arma::Cube<ValueType>& result) const {
        arma::Cube<ValueType> scalarFunctionValues;
        ScalarSpaceMapping<ValueType>::
                evaluateShapeFunctions(basisData, geomData, scalarFunctionValues);
        const arma::Mat<CoordinateType>& n = geomData.normals;

        const size_t dimWorld = 3;
        const size_t functionCount = scalarFunctionValues.n_cols;
        const size_t pointCount = scalarFunctionValues.n_slices;
        result.set_size(dimWorld, functionCount, pointCount);
        for (size_t p = 0; p < pointCount; ++p)
            for (size_t f = 0; f < functionCount; ++f)
                for (size_t d = 0; d < dimWorld; ++d)
                    result(d, f, p) = scalarFunctionValues(0, f, p) * n(d, p);
    }

    virtual void evaluateImplReal(const BasisData<CoordinateType>& basisData,
                                  const GeometricalData<CoordinateType>& geomData,
                                  arma::Cube<CoordinateType>& result) const {
        evaluateImpl(basisData, geomData, result);
    }

    virtual void evaluateImplComplex(const BasisData<ComplexType>& basisData,
                                     const GeometricalData<CoordinateType>& geomData,
                                     arma::Cube<ComplexType>& result) const {
        evaluateImpl(basisData, geomData, result);
    }

};

} // namespace Fiber

#endif
