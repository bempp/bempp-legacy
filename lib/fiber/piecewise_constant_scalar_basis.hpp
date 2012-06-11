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

#ifndef fiber_piecewise_constant_scalar_basis_hpp
#define fiber_piecewise_constant_scalar_basis_hpp

#include "basis.hpp"
#include "basis_data.hpp"

#include <algorithm>

namespace Fiber
{

template <typename ValueType>
class PiecewiseConstantScalarBasis : public Basis<ValueType>
{
public:
    typedef typename Basis<ValueType>::CoordinateType CoordinateType;

    virtual size_t size() const {
        return 1;
    }

    virtual size_t order() const {
        return 0;
    }

    virtual void evaluate(size_t what,
                          const arma::Mat<CoordinateType>& points,
                          LocalDofIndex localDofIndex,
                          BasisData<ValueType>& data) const {
        if (localDofIndex != ALL_DOFS && localDofIndex != 0)
            throw std::invalid_argument("PiecewiseConstantScalarBasis::evaluate(): "
                                        "Invalid localDofIndex");
        // Since there is only one basis function, there is no difference
        // between calculating all basis functions and just one.

        const int componentCount = 1;
        const int functionCount = 1;
        const int pointCount = points.n_cols;
        if (what & VALUES)
        {
            data.values.set_size(componentCount, functionCount, pointCount);
            data.values.fill(1.);
        }
        if (what & DERIVATIVES)
        {
            const int coordCount = points.n_rows;
            data.derivatives.set_size(componentCount, coordCount,
                                 functionCount, pointCount);
            std::fill(data.derivatives.begin(), data.derivatives.end(), 0.);
        }
    }
};

} // namespace Fiber

#endif
