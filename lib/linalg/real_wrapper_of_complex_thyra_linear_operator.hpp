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

#include "config_trilinos.hpp"

#ifndef bempp_real_wrapper_of_complex_thyra_linear_operator_hpp
#define bempp_real_wrapper_of_complex_thyra_linear_operator_hpp

#include "../common/common.hpp"


#ifdef WITH_TRILINOS
#include <Thyra_LinearOpDefaultBase_decl.hpp>

namespace Bempp
{

template <typename ValueType>
class RealWrapperOfComplexThyraLinearOperator
        : public Thyra::LinearOpDefaultBase<ValueType>
{
public:
    typedef std::complex<ValueType> ComplexValueType;
    typedef Thyra::LinearOpBase<std::complex<ValueType> > ComplexLinearOp;

    RealWrapperOfComplexThyraLinearOperator(
            const Teuchos::RCP<const ComplexLinearOp>& complexOperator);

    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > domain() const;
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > range() const;

protected:
    virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const;

    virtual void applyImpl(
            const Thyra::EOpTransp M_trans,
            const Thyra::MultiVectorBase<ValueType> &X_in,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
            const ValueType alpha,
            const ValueType beta) const;

private:
    Teuchos::RCP<const ComplexLinearOp> m_complexOperator;
    Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > m_domainSpace;
    Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> > m_rangeSpace;
};

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
