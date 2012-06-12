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

#ifndef fiber_scalar_space_mapping_hpp
#define fiber_scalar_space_mapping_hpp

#include "../common/common.hpp"

#include <cassert>
#include <stdexcept>
#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

namespace Fiber
{

/** \brief Mapping of functions defined on reference elements to ones defined
  on physical elements. */
template <typename ValueType>
class ScalarSpaceMapping
{
public:
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    static void addShapeFunctionDependencies(size_t& basisDeps, size_t& geomDeps) {
        basisDeps |= VALUES;
    }

    static void addSurfaceCurlDependencies(size_t& basisDeps, size_t& geomDeps) {
        basisDeps |= DERIVATIVES;
        geomDeps |= NORMALS | JACOBIAN_INVERSES_TRANSPOSED;
    }

    static void evaluateShapeFunctions(const BasisData<ValueType>& basisData,
                                       const GeometricalData<CoordinateType>& geomData,
                                       arma::Cube<ValueType>& result) {
        result = basisData.values;
    }

    /**
      \param[out] result
        dimensions: (worldDim, functionCount, pointCount)
    */
    static void evaluateSurfaceCurls3d(const BasisData<ValueType>& basisData,
                                       const GeometricalData<CoordinateType>& geomData,
                                       arma::Cube<ValueType>& result) {
        const arma::Mat<CoordinateType>& n = geomData.normals;
        // jt(i, j): dx_j/dq_i
        const arma::Cube<CoordinateType>& jit = geomData.jacobianInversesTransposed;
        const Array4d<ValueType>& d = basisData.derivatives;

        assert(d.extent(0) == 1); // scalar functions
        const int worldDim = 3;
        assert(static_cast<int>(d.extent(1)) + 1 == worldDim);
        const size_t functionCount = d.extent(2);
        const size_t pointCount = d.extent(3);

        result.set_size(worldDim, functionCount, pointCount);

        for (size_t p = 0; p < pointCount; ++p)
            for (size_t f = 0; f < functionCount; ++f) {
                // vec := gradient of the basis function extended outside
                // the surface so that its normal derivative on the surf. is zero
                arma::Col<ValueType> vec(3);
                vec(0) = d(0, 0, f, p) * jit(0, 0, p) + d(0, 1, f, p) * jit(0, 1, p);
                vec(1) = d(0, 0, f, p) * jit(1, 0, p) + d(0, 1, f, p) * jit(1, 1, p);
                vec(2) = d(0, 0, f, p) * jit(2, 0, p) + d(0, 1, f, p) * jit(2, 1, p);

                // result := n \times vec
                result(0, f, p) = n(1, p) * vec(2) - n(2, p) * vec(1);
                result(1, f, p) = n(2, p) * vec(0) - n(0, p) * vec(2);
                result(2, f, p) = n(0, p) * vec(1) - n(1, p) * vec(0);
            }
    }

};

} // namespace Fiber

#endif // SCALAR_SPACE_MAPPING_HPP
