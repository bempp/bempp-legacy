// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_scalar_function_value_times_normal_functor_hpp
#define fiber_scalar_function_value_times_normal_functor_hpp

#include "../common/common.hpp"

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "basis_transformation_functor_wrappers.hpp"

namespace Fiber
{

template <typename CoordinateType_>
class ScalarFunctionValueTimesNormalElementaryFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int argumentDimension() const { return 1; }
    int resultDimension() const { return 3; }

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        basisDeps |= VALUES;
        geomDeps |= NORMALS;
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            _1dSliceOf3dArray<ValueType>& result) const {
        assert(basisData.componentCount() == 1);
        const int dimWorld = geomData.dimWorld();
        for (int dim = 0; dim < dimWorld; ++dim)
            result(dim) = basisData.values(0) * geomData.normal(dim);
    }
};

// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
/** \ingroup functors
 *  \brief Functor calculating the value of a scalar shape function multiplied
 *  by the unit vector normal to the surface. */
template <typename CoordinateType_>
class ScalarFunctionValueTimesNormalFunctor :
        public ElementaryBasisTransformationFunctorWrapper<
        ScalarFunctionValueTimesNormalElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
