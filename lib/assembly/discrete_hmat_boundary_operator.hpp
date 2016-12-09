// Copyright (C) 2011-2014 by the BEM++ Authors
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

#ifndef bempp_discrete_hmat_boundary_operator_hpp
#define bempp_discrete_hmat_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"
#include "../common/shared_ptr.hpp"
#include "discrete_boundary_operator.hpp"
#include "../hmat/hmatrix.hpp"

namespace Bempp {

template <typename ValueType>
class DiscreteHMatBoundaryOperator
    : public DiscreteBoundaryOperator<ValueType> {
public:
  DiscreteHMatBoundaryOperator(
      const shared_ptr<hmat::DefaultHMatrixType<ValueType>> &hMatrix);

  unsigned int rowCount() const override;

  unsigned int columnCount() const override;

  shared_ptr<const hmat::DefaultHMatrixType<ValueType>> hMatrix() const;

  void addBlock(const std::vector<int> &rows, const std::vector<int> &cols,
                const ValueType alpha, Matrix<ValueType> &block) const override;

private:
  void applyBuiltInImpl(const TranspositionMode trans,
                        const Eigen::Ref<Vector<ValueType>> &x_in,
                        Eigen::Ref<Vector<ValueType>> y_inout,
                        const ValueType alpha,
                        const ValueType beta) const override;

  shared_ptr<hmat::DefaultHMatrixType<ValueType>> m_hMatrix;
};

template <typename ValueType>
shared_ptr<const hmat::DefaultHMatrixType<ValueType>>
castToHMatrix(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op);
}

#endif
