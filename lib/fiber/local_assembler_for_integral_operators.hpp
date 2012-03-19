// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_local_assembler_for_integral_operators_hpp
#define fiber_local_assembler_for_integral_operators_hpp

#include "array_2d.hpp"
#include "types.hpp"
#include <vector>
#include <stdexcept>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType> class TestKernelTrialIntegrator;

// Things to keep in mind:
//
// Additional quadrature order specific to the Helmholtz equation (dependence
// on k*h) will need to be handled via a kernel-specific code (there will be a
// particular quadrature selector for Helmholtz kernels, probably derived from
// the standard Sauter-Schwab quadrature selector).
//
// The quadrature manager needs to know the polynomial order of an element
// and the increment in the integration order resulting from mesh curvature.

/** \brief Local assembler for integral operators.

  An integration manager provides methods that choose appropriate quadrature
  rules (weights and points) for integrals occurring in boundary-element
  matrices of integral operators. Factors influencing the choice of quadrature
  rule include the regularity properties of a kernel, maximum polynomial order
  of basis functions on an element, element volume and (for double integrals)
  element-element distance.
 */
template <typename ValueType>
class LocalAssemblerForIntegralOperators
{
public:
    virtual ~LocalAssemblerForIntegralOperators() {}

//    /** \brief Number of dimensions over which integration takes place. */
//    int elementDimension() const {
//        asImp().elementDimension();
//    }    

    /** \brief Assemble local weak forms.

    In this overload, a "column" of local weak forms is assembled. More
    specifically, on exit \p result is a vector of local weak forms corresponding
    to the following pairs of elements:

    - if \p callVariant is \p TEST_TRIAL, all pairs (\p elementA, \p elementB)
    for \p elementA in \p elementsA;

    - if \p callVariant is \p TRIAL_TEST, all pairs (\p elementB, \p elementA)
    for \p elementA in \p elementsA.

    Unless \p localDofIndexB is set to \p ALL_DOFS, only entries corresponding
    to the (\p localDofIndexB)th local DOF on \p elementB are calculated. */
    virtual void evaluateLocalWeakForms(
            CallVariant callVariant,
            const std::vector<int>& elementIndicesA,
            int elementIndexB,
            LocalDofIndex localDofIndexB,
            std::vector<arma::Mat<ValueType> >& result) = 0;

    /** \brief Assemble local weak forms.

    This overload constructs and assigns to the output parameter \p result the
    2D array of local weak forms corresponding to all pairs (\p testElement, \p
    trialElement) with \p testElement in \p testElementIndices and \p
    trialElement in \p trialElementIndices.

    This function should be used primarily for small blocks of elements lying
    close to each other. */
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& testElementIndices,
            const std::vector<int>& trialElementIndices,
            Fiber::Array2D<arma::Mat<ValueType> >& result) = 0;

    // TODO: evaluateLocalOperator or something similar
};

} // namespace Fiber

#endif
