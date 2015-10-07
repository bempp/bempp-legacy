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

#ifndef fiber_raviart_thomas_0_shapeset_barycentric_hpp
#define fiber_raviart_thomas_0_shapeset_barycentric_hpp

#include "basis.hpp"
#include "basis_data.hpp"
#include "raviart_thomas_0_shapeset.hpp"

#include <math.h>

namespace Fiber {

template <typename ValueType>
class RaviartThomas0ShapesetBarycentric : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;
  enum BasisType { TYPE1, TYPE2 };

public:
  RaviartThomas0ShapesetBarycentric(BasisType type) : m_type(type) {}

  virtual int size() const { return 3; }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {

    BasisData<ValueType> temp;

    raviartBasis.evaluate(what, points, ALL_DOFS, temp);

    ValueType coeffs[2][3][3] = {
        {{0,0,0}, {0, 0, 0}, {0, 0, sqrt(58)/12}},
        {{0,0,0}, {0, 0, 0}, {0, 0, sqrt(58)/12}}};

    int type_index;
    if (m_type == TYPE1)
      type_index = 0;
    else
      type_index = 1;

    if (localDofIndex != ALL_DOFS) {

      if (what & VALUES) {
        data.values.set_size(1, 1, temp.values.extent(2));
        for (int i = 0; i < data.values.extent(2); ++i) {
          data.values(0, 0, i) = 0;
          for (int j = 0; j < 3; ++j)
            data.values(0, 0, i) +=
                coeffs[type_index][localDofIndex][j] * temp.values(0, j, i);
        }
      }
      if (what & DERIVATIVES) {
        data.derivatives.set_size(1, temp.derivatives.extent(1), 1,
                                  temp.derivatives.extent(3));
        for (int i = 0; i < data.derivatives.extent(1); ++i)
          for (int j = 0; j < data.derivatives.extent(3); ++j) {
            data.derivatives(0, i, 0, j) = 0;
            for (int k = 0; k < 3; ++k)
              data.derivatives(0, i, 0, j) +=
                  coeffs[type_index][localDofIndex][k] *
                  temp.derivatives(0, i, k, j);
          }
      }

    } else {
      if (what & VALUES) {
        data.values.set_size(1, 3, temp.values.extent(2));
        for (int dofIndex = 0; dofIndex < 3; ++dofIndex) {
          for (int i = 0; i < data.values.extent(2); ++i) {
            data.values(0, dofIndex, i) = 0;
            for (int j = 0; j < 3; ++j)
              data.values(0, dofIndex, i) +=
                  coeffs[type_index][dofIndex][j] * temp.values(0, j, i);
          }
        }
      }

      if (what & DERIVATIVES) {
        data.derivatives.set_size(1, temp.derivatives.extent(1), 3,
                                  temp.derivatives.extent(3));
        for (int dofIndex = 0; dofIndex < 3; ++dofIndex)
          for (int i = 0; i < data.derivatives.extent(1); ++i)
            for (int j = 0; j < data.derivatives.extent(3); ++j) {
              data.derivatives(0, i, dofIndex, j) = 0;
              for (int k = 0; k < 3; ++k)
                data.derivatives(0, i, dofIndex, j) +=
                    coeffs[type_index][dofIndex][k] *
                    temp.derivatives(0, i, k, j);
            }
      }
    }
  }

  virtual std::pair<const char *, int> clCodeString(bool isTestBasis) const {
    throw std::runtime_error(
        "RaviartThomas0BasisBarycentric::clCodeString():"
        "OpenCL not supported for this basis type.");
  }

private:
  Fiber::RaviartThomas0Shapeset<3, ValueType> raviartBasis;
  BasisType m_type;
};

} // namespace Fiber

#endif

