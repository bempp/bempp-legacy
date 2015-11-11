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
//  RaviartThomas0ShapesetBarycentric(BasisType type, int testInt) : m_type(type) {//std::cout << testInt << std::endl;}
  RaviartThomas0ShapesetBarycentric(){}
  RaviartThomas0ShapesetBarycentric(BasisType type) : m_type(type){
    Matrix<double> coeffs;
    coeffs.conservativeResize(3,3);

    if(type==TYPE1){
        coeffs(0,0)=1./3.;     coeffs(0,1)=0.;        coeffs(0,2)=-1./6.;
        coeffs(1,0)=-1./3.;    coeffs(1,1)=1./2.;     coeffs(1,2)=0.;
        coeffs(2,0)=0.;        coeffs(2,1)=0.;        coeffs(2,2)=1./6.;

    } else {
        coeffs(0,0)=1./2.;     coeffs(0,1)=-1./3.;    coeffs(0,2)=0.;
        coeffs(1,0)=0.;        coeffs(1,1)=1./3.;     coeffs(1,2)=-1./6.;
        coeffs(2,0)=0.;        coeffs(2,1)=0.;        coeffs(2,2)=1./6.;

    }
    m_coeffs = coeffs.cast<ValueType>();
  }

  virtual int size() const { return 3; }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {

    BasisData<ValueType> temp;
    raviartBasis.evaluate(what, points, ALL_DOFS, temp);

    if (localDofIndex != ALL_DOFS) {

      if (what & VALUES) {
        data.values.set_size(1, 1, temp.values.extent(2));
        for (int i = 0; i < data.values.extent(2); ++i) {
          for (int k=0; k < data.values.extent(0); ++k){
            data.values(k, 0, i) = 0;
            for (int j = 0; j < 3; ++j)
              data.values(k, 0, i) +=
                  m_coeffs(localDofIndex,j) * temp.values(k, j, i);
          }
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
                  m_coeffs(localDofIndex,k) *
                  temp.derivatives(0, i, k, j);
          }
      }

    } else {
      if (what & VALUES) {
        data.values.set_size(temp.values.extent(0), 3, temp.values.extent(2));
        for (int dofIndex = 0; dofIndex < 3; ++dofIndex) {
          for (int i = 0; i < data.values.extent(2); ++i) {
            for (int k=0; k < data.values.extent(0); ++k){
              data.values(k, dofIndex, i) = 0;
              for (int j = 0; j < 3; ++j)
                data.values(k, dofIndex, i) +=
                    m_coeffs(dofIndex,j) * temp.values(k, j, i);
            }
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
                    m_coeffs(dofIndex,k) *
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
  mutable BasisType m_type;
  mutable int m_sonIndex;
  Matrix<ValueType> m_coeffs;
};

} // namespace Fiber

#endif

