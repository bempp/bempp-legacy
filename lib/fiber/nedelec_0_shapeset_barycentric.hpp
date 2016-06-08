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

#ifndef fiber_nedelec_0_shapeset_barycentric_hpp
#define fiber_nedelec_0_shapeset_barycentric_hpp

#include "basis.hpp"
#include "basis_data.hpp"
#include "raviart_thomas_0_shapeset_barycentric.hpp"

namespace Fiber {

template <typename ValueType>
class Nedelec0ShapesetBarycentric : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;
  enum BasisType { TYPE1, TYPE2 };

public:
  Nedelec0ShapesetBarycentric(BasisType type) : m_type(type){}

  virtual int size() const { return 3; }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {

    BasisData<ValueType> temp;
    if (what & VALUES) {
        if (m_type == TYPE1) raviartBasis1.evaluate(what, points, localDofIndex, temp);
        else                 raviartBasis2.evaluate(what, points, localDofIndex, temp);

        data.values.set_size(temp.values.extent(0),temp.values.extent(1),temp.values.extent(2));
        for (int i=0; i!=temp.values.extent(1); ++i)
          for (int j=0; j!=temp.values.extent(2); ++j) {
            data.values(0,i,j) = -temp.values(1,i,j);
            data.values(1,i,j) = temp.values(0,i,j);
          }
    }
    if (what & DERIVATIVES) {
        if (m_type == TYPE1) raviartBasis1.evaluate(what, points, localDofIndex, temp);
        else                 raviartBasis2.evaluate(what, points, localDofIndex, temp);

        data.derivatives.set_size(temp.derivatives.extent(0),temp.derivatives.extent(1),temp.derivatives.extent(2),temp.derivatives.extent(3));
        for (int i=0; i!=temp.derivatives.extent(1); ++i)
          for (int j=0; j!=temp.derivatives.extent(2); ++j)
            for (int k=0; k!=temp.derivatives.extent(3); ++k) {
              data.derivatives(0,i,j,k) = -temp.derivatives(1,i,j,k);
              data.derivatives(1,i,j,k) = temp.derivatives(0,i,j,k);
            }
    }

  }

  virtual std::pair<const char *, int> clCodeString(bool isTestBasis) const {
    throw std::runtime_error(
        "Nedelec0BasisBarycentric::clCodeString():"
        "OpenCL not supported for this basis type.");
  }

private:
  typedef Fiber::RaviartThomas0ShapesetBarycentric<ValueType> RTShapeset;
  RTShapeset raviartBasis1 = RTShapeset(RTShapeset::TYPE1);
  RTShapeset raviartBasis2 = RTShapeset(RTShapeset::TYPE2);
  mutable BasisType m_type;
};

} // namespace Fiber

#endif
