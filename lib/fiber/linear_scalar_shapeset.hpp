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

#ifndef fiber_linear_scalar_shapeset_hpp
#define fiber_linear_scalar_shapeset_hpp

#include "../common/common.hpp"

#include "basis.hpp"

#include "basis_data.hpp"
#include "dune_basis_helper.hpp"
#include "CL/piecewise_linear_continuous_scalar_basis.cl.str"

#include <dune/localfunctions/lagrange/p1/p1localbasis.hh>
#include <dune/localfunctions/lagrange/q1/q1localbasis.hh>

namespace Fiber {

template <int elementVertexCount, typename CoordinateType, typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits {};

// Line
template <typename CoordinateType, typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<2, CoordinateType,
                                                  ValueType> {
public:
  typedef Dune::Q1LocalBasis<CoordinateType, ValueType, 1> DuneBasis;
};

// Triangle
template <typename CoordinateType, typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<3, CoordinateType,
                                                  ValueType> {
public:
  typedef Dune::P1LocalBasis<CoordinateType, ValueType, 2> DuneBasis;
};

// Quadrilateral
template <typename CoordinateType, typename ValueType>
struct PiecewiseLinearContinuousScalarBasisTraits<4, CoordinateType,
                                                  ValueType> {
public:
  typedef Dune::Q1LocalBasis<CoordinateType, ValueType, 2> DuneBasis;
};

/** \brief Shapeset composed of the linear nodal functions. */
template <int elementVertexCount, typename ValueType>
class LinearScalarShapeset : public Basis<ValueType> {
public:
  typedef typename Basis<ValueType>::CoordinateType CoordinateType;

private:
  typedef typename PiecewiseLinearContinuousScalarBasisTraits<
      elementVertexCount, CoordinateType, ValueType>::DuneBasis DuneBasis;

public:
  virtual int size() const {
    DuneBasis basis;
    return basis.size();
  }

  virtual int order() const { return 1; }

  virtual void evaluate(size_t what, const Matrix<CoordinateType> &points,
                        LocalDofIndex localDofIndex,
                        BasisData<ValueType> &data) const {
    if (localDofIndex != ALL_DOFS &&
        (localDofIndex < 0 || size() <= localDofIndex))
      throw std::invalid_argument("LinearScalarShapeset::"
                                  "evaluate(): Invalid localDofIndex");

    if (what & VALUES)
      evaluateShapeFunctionsWithDune<CoordinateType, ValueType, DuneBasis>(
          points, localDofIndex, data.values);
    if (what & DERIVATIVES)
      evaluateShapeFunctionDerivativesWithDune<CoordinateType, ValueType,
                                               DuneBasis>(points, localDofIndex,
                                                          data.derivatives);
  }

  virtual std::pair<const char *, int> clCodeString(bool isTestShapeset) const {
    static struct StrBuf {
      std::pair<char *, int> str;
      bool need_setup;
    } test = {std::pair<char *, int>(NULL, 0), true},
      trial = {std::pair<char *, int>(NULL, 0), true};

    StrBuf &buf = (isTestShapeset ? test : trial);
    const char *modifier = (isTestShapeset ? "A" : "B");

    if (buf.need_setup) {
      std::string funcName("clBasisf");
      std::string code(piecewise_linear_continuous_scalar_basis_cl,
                       piecewise_linear_continuous_scalar_basis_cl_len);
      size_t len = funcName.size();
      size_t n = code.find(funcName);
      while (n != std::string::npos) {
        code.insert(n + len, modifier);
        n = code.find(funcName, n + len);
      }
      buf.str.second = code.size();
      buf.str.first = new char[buf.str.second + 1];
      strcpy(buf.str.first, code.c_str());
      buf.need_setup = false;
    }
    return buf.str;
  }
};

} // namespace Fiber

#endif
