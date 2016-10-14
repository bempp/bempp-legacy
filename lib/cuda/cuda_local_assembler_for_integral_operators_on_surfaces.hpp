// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_cuda_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_cuda_local_assembler_for_integral_operators_on_surfaces_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename ResultType> class ScalarTraits;
/** \endcond */


template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class CudaLocalAssemblerForIntegralOperatorsOnSurfaces
    : public LocalAssemblerForIntegralOperators<ResultType> {

public:

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  CudaLocalAssemblerForIntegralOperatorsOnSurfaces();

  virtual ~CudaLocalAssemblerForIntegralOperatorsOnSurfaces();

  virtual void evaluateLocalWeakForms(CallVariant callVariant,
                                      const std::vector<int> &elementIndicesA,
                                      int elementIndexB,
                                      LocalDofIndex localDofIndexB,
                                      std::vector<Matrix<ResultType>> &result,
                                      CoordinateType nominalDistance = -1.);

};
} // namespace Fiber

#endif
