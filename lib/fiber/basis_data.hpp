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

#ifndef fiber_basis_data_hpp
#define fiber_basis_data_hpp

#include <armadillo>
#include "../common/multidimensional_arrays.hpp"

namespace Fiber
{

enum BasisDataType
{
    VALUES = 0x0001,
    DERIVATIVES = 0x0002
};

template <typename ValueType> class ConstBasisDataSlice;

template <typename ValueType>
struct BasisData
{
    // values(i,k,l) = (f_k)_i(x_l) ->
    // ith component of kth basis function at lth point
    arma::Cube<ValueType> values;
    // derivatives(i,j,k,l) = (d_j (f_k)_i)(x_l) ->
    // derivative in direction j of ith component of kth basis function at lth point
    _4dArray<ValueType> derivatives;

    int componentCount() const {
        return std::max<int>(values.n_rows, derivatives.extent(0));
    }

    int functionCount() const {
        return std::max<int>(values.n_cols, derivatives.extent(2));
    }

    int pointCount() const {
        return std::max<int>(values.n_slices, derivatives.extent(3));
    }

    ConstBasisDataSlice<ValueType> const_slice(int function, int point) const {
        return ConstBasisDataSlice<ValueType>(*this, function, point);
    }
};

template <typename ValueType>
class ConstBasisDataSlice
{
public:
    ConstBasisDataSlice(const BasisData<ValueType>& basisData,
                        int function, int point) :
        m_basisData(basisData), m_function(function), m_point(point) {}

    ValueType values(int dim) const {
        return m_basisData.values(dim, m_function, m_point);
    }

    ValueType derivatives(int dim, int direction) const {
        return m_basisData.derivatives(dim, direction, m_function, m_point);
    }

    int componentCount() const {
        return m_basisData.componentCount();
    }

private:
    const BasisData<ValueType>& m_basisData;
    int m_function, m_point;
};

} // namespace Fiber

#endif // BASIS_DATA_TYPES_HPP
