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

#ifndef bempp_elementary_integral_operator_hpp
#define bempp_elementary_integral_operator_hpp

#include "linear_operator.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/types.hpp"
#include "../fiber/kernel.hpp"
#include "../fiber/types.hpp"

#include <vector>
#include <armadillo>

namespace Fiber
{

template <typename ValueType> class Expression;
template <typename ValueType> class LocalAssemblerForIntegralOperators;

} // namespace Fiber

namespace Bempp
{

template <int codim> class EntityPointer;

template <typename ValueType> class WeakFormAcaAssemblyHelper;

template <typename ValueType>
class ElementaryIntegralOperator : public LinearOperator<ValueType>
{
public:
    typedef typename LinearOperator<ValueType>::LocalAssemblerFactory
    LocalAssemblerFactory;
    typedef typename Fiber::LocalAssemblerForIntegralOperators<ValueType>
    LocalAssembler;

    virtual int trialComponentCount() const {
        return kernel().domainDimension();
    }

    virtual int testComponentCount() const {
        return kernel().codomainDimension();
    }

    virtual bool isRegular() const = 0;

    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperator(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

private:
    virtual const Fiber::Kernel<ValueType>& kernel() const = 0;
    virtual const Fiber::Expression<ValueType>& testExpression() const = 0;
    virtual const Fiber::Expression<ValueType>& trialExpression() const = 0;

    /** @}
        \name Operator assembly
        @{ */
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInDenseMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInAcaMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >

    /** @}
        \name Weak form assembly
        @{ */
    assembleWeakFormInDenseMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInAcaMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
    /** @} */
};

} // namespace Bempp

#endif
