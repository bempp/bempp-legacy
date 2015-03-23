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

#ifndef fiber_local_assembler_for_local_operators_hpp
#define fiber_local_assembler_for_local_operators_hpp

#include "../common/common.hpp"

#include "_2d_array.hpp"
#include "scalar_traits.hpp"
#include "types.hpp"

#include <vector>

namespace Fiber {

/** \brief Abstract interface of a local assembler for local operators.

  Let \f$L\f$ be a local boundary operator. A local assembler associated with
  this operator constructs matrices \f$A_e\f$, where \f$e\f$ is an element
  number, with entries \f$(A_{e})_{ij} = (u_{ei}, A v_{ej}\f$) where
  \f$(\cdot, \cdot)\f$ is the standard scalar product, \f$(u_{ei})_{i=1}^m\f$
  are the element-level test functions defined on element \f$e\f$ and
  \f$(v_{ej})_{j=1}^n\f$ are the element-level trial functions defined on that
  element.

  The local assembler is responsible for choosing an appropriate way of
  evaluating
  the necessary integrals. */
template <typename ResultType> class LocalAssemblerForLocalOperators {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  virtual ~LocalAssemblerForLocalOperators() {}

  /** \brief Assemble local weak forms.

    \param[in]  elementIndices Vector of element indices.
    \param[out] result         Vector of weak forms of the operator on
                               element pairs
                               (\p element(\p i), \p element(\p i))
                               for \p i in \p elementIndices. */
  virtual void
  evaluateLocalWeakForms(const std::vector<int> &elementIndices,
                         std::vector<Matrix<ResultType>> &result) = 0;
};

} // namespace Fiber

#endif
