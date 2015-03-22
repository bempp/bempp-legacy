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

#ifndef fiber_constant_scalar_shapeset_hpp
#define fiber_constant_scalar_shapeset_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "basis.hpp"
#include "basis_data.hpp"

#include <algorithm>

namespace Fiber {

/** \brief A shapeset containing only one function: the constant function. */
template <typename ValueType>
class ConstantScalarShapeset : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;

  virtual int size() const { return 1; }

  virtual int order() const { return 0; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {
    if (localDofIndex != ALL_DOFS && localDofIndex != 0)
      throw std::invalid_argument("ConstantScalarShapeset::evaluate(): "
                                  "Invalid localDofIndex");
    // Since there is only one shape function, there is no difference
    // between calculating all shape functions and just one.

    const int componentCount = 1;
    const int functionCount = 1;
    const int pointCount = points.n_cols;
    if (what & VALUES) {
      data.values.set_size(componentCount, functionCount, pointCount);
      std::fill(data.values.begin(), data.values.end(), 1.);
    }
    if (what & DERIVATIVES) {
      const int coordCount = points.n_rows;
      data.derivatives.set_size(componentCount, coordCount, functionCount,
                                pointCount);
      std::fill(data.derivatives.begin(), data.derivatives.end(), 0.);
    }
  }
};

} // namespace Fiber

#endif
