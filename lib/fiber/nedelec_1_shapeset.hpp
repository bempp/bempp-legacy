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

#ifndef fiber_nedelec_1_shapeset_hpp
#define fiber_nedelec_1_shapeset_hpp

#include "basis.hpp"
#include "basis_data.hpp"
#include "raviart_thomas_0_shapeset.hpp"

namespace Fiber {

template <int elementVertexCount, typename CoordinateType, typename ValueType>
struct Nedelec1BasisTraits {};

// Triangle
template <typename CoordinateType, typename ValueType>
struct Nedelec1BasisTraits<3, CoordinateType, ValueType> {
public:
  typedef Dune::RT02DLocalBasis<CoordinateType, ValueType> DuneBasis;
};

/** \brief Shapeset composed of the lowest-order Raviart-Thomas functions. */
template <int elementVertexCount, typename ValueType>
class Nedelec1Shapeset : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;

private:
  typedef typename Nedelec1BasisTraits<elementVertexCount, CoordinateType,
                                             ValueType>::DuneBasis DuneBasis;

public:
  virtual int size() const {
    DuneBasis basis;
    return basis.size();
  }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {


    if (what & VALUES) {
        BasisData<ValueType> temp;
        raviartBasis.evaluate(what, points, ALL_DOFS, temp);
        data.values.set_size(temp.values.extent(0),temp.values.extent(1),temp.values.extent(2));
        for (int dofIndex = 0; dofIndex != temp.values.extent(1); ++dofIndex)
            for (int i=0; i != temp.values.extent(2); ++i){
                data.values(0,dofIndex,i) = temp.values(1,dofIndex,i);
                data.values(1,dofIndex,i) = -temp.values(0,dofIndex,i);
            }
    }
    if (what & DERIVATIVES) {
        raviartBasis.evaluate(what, points, ALL_DOFS, data);
    }

  }
private:
    Fiber::RaviartThomas0Shapeset<3, ValueType> raviartBasis;

};

} // namespace Fiber

#endif
