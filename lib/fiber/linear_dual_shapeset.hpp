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

#ifndef fiber_linear_dual_shapeset_hpp
#define fiber_linear_dual_shapeset_hpp

#include "basis.hpp"
#include "basis_data.hpp"
#include "linear_scalar_shapeset.hpp"

namespace Fiber {

template <typename ValueType>
class LinearDualShapeset : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;

public:
  LinearDualShapeset() {}
  LinearDualShapeset(Matrix<ValueType> coeffs) : m_coeffs(coeffs) {}

  virtual int size() const { return m_coeffs.cols(); }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {

    BasisData<ValueType> temp;
    linearBasis.evaluate(what, points, ALL_DOFS, temp);

    if (localDofIndex != ALL_DOFS) {
      if (what & VALUES) {
        data.values.set_size(1, 1, temp.values.extent(2));
        for (int i = 0; i < data.values.extent(2); ++i) {
          data.values(0, 0, i) = 0;
          for (int j = 0; j < 3; ++j)
            data.values(0, 0, i) += m_coeffs(j, localDofIndex) * temp.values(0, j, i);
          }
      }
      if (what & DERIVATIVES) {
        data.derivatives.set_size(1,
                                  temp.derivatives.extent(1), 1,
                                  temp.derivatives.extent(3));
        for (int i = 0; i < data.derivatives.extent(1); ++i)
          for (int j = 0; j < data.derivatives.extent(3); ++j) {
            data.derivatives(0, i, 0, j) = 0;
            for (int k = 0; k < 3; ++k)
              data.derivatives(0, i, 0, j) += m_coeffs(k, localDofIndex) * temp.derivatives(0, i, k, j);
            }
      }
    } else {
      if (what & VALUES) {
        data.values.set_size(1, m_coeffs.cols(), temp.values.extent(2));
        for (int dofIndex = 0; dofIndex < m_coeffs.cols(); ++dofIndex)
          for (int i = 0; i < data.values.extent(2); ++i) {
            data.values(0, dofIndex, i) = 0;
            for (int j = 0; j < 3; ++j)
              data.values(0, dofIndex, i) += m_coeffs(j, dofIndex) * temp.values(0, j, i);
          }
      }
      if (what & DERIVATIVES) {
        data.derivatives.set_size(1,
                                  temp.derivatives.extent(1), m_coeffs.cols(),
                                  temp.derivatives.extent(3));
        for (int dofIndex = 0; dofIndex < m_coeffs.cols(); ++dofIndex)
          for (int i = 0; i < data.derivatives.extent(1); ++i)
            for (int j = 0; j < data.derivatives.extent(3); ++j) {
              data.derivatives(0, i, dofIndex, j) = 0;
              for (int k = 0; k < 3; ++k)
                data.derivatives(0, i, dofIndex, j) += m_coeffs(k, dofIndex) * temp.derivatives(0, i, k, j);
            }
      }
    }
  }

  virtual std::pair<const char *, int> clCodeString(bool isTestBasis) const {
    throw std::runtime_error("LinearDualBasis::clCodeString():"
                             "OpenCL not supported for this basis type.");
  }

private:
  Fiber::LinearScalarShapeset<3, ValueType> linearBasis;
  mutable Matrix<ValueType> m_coeffs;
};

} // namespace Fiber

#endif
