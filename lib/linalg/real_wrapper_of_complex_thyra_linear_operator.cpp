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

#include "real_wrapper_of_complex_thyra_linear_operator.hpp"
#ifdef WITH_TRILINOS

#include "../fiber/explicit_instantiation.hpp"

#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_DefaultSpmdVectorSpace.hpp>
#include <Thyra_VectorSpaceBase.hpp>

namespace Bempp
{

template <typename ValueType>
RealWrapperOfComplexThyraLinearOperator<ValueType>::
RealWrapperOfComplexThyraLinearOperator(
        const Teuchos::RCP<const ComplexLinearOp>& complexOperator) :
    m_complexOperator(complexOperator),
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(
                      2 * complexOperator->domain()->dim())),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(
                     2 * complexOperator->range()->dim()))
{
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
RealWrapperOfComplexThyraLinearOperator<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
RealWrapperOfComplexThyraLinearOperator<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool RealWrapperOfComplexThyraLinearOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    return m_complexOperator->opSupported(M_trans);
}

template <typename ValueType>
void RealWrapperOfComplexThyraLinearOperator<ValueType>::applyImpl(
        const Thyra::EOpTransp M_trans,
        const Thyra::MultiVectorBase<ValueType> &X_in,
        const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    const Teuchos::Ordinal rowCount_X_in = X_in.range()->dim();
    assert(rowCount_X_in % 2 == 0); // each complex number has two parts
    const Teuchos::Ordinal colCount_X_in = X_in.domain()->dim();

    Teuchos::ArrayRCP<const ComplexValueType> complexArray_X_in;
    RTOpPack::ConstSubMultiVectorView<ComplexValueType> complexView_X_in;
    Teuchos::RCP<const Thyra::MultiVectorBase<ComplexValueType> > complex_X_in;

    Thyra::ConstDetachedMultiVectorView<ValueType> view_X_in(
                Teuchos::rcpFromRef(X_in));
    if (view_X_in.leadingDim() == rowCount_X_in) {
        // contiguous vector
        const ComplexValueType* complexData_X_in =
                reinterpret_cast<const ComplexValueType*>(view_X_in.values());
        complexArray_X_in = Teuchos::arcp(complexData_X_in,
                                          0, // lowerOffset
                                          (rowCount_X_in / 2)* colCount_X_in,
                                          false); // doesn't own
        complexView_X_in.initialize(0, // globalOffset
                                    rowCount_X_in / 2,
                                    0, // colOffset
                                    colCount_X_in,
                                    complexArray_X_in,
                                    rowCount_X_in / 2); // leadingDim
        complex_X_in = Thyra::createMembersView(m_complexOperator->domain(),
                                                complexView_X_in);
    } else
        throw std::runtime_error("RealWrapperOfComplexThyraLinearOperator::"
                                 "applyImpl(): discontiguous multivectors are not "
                                 "supported yet");

    const Teuchos::Ordinal rowCount_Y_inout = Y_inout->range()->dim();
    assert(rowCount_Y_inout % 2 == 0);
    const Teuchos::Ordinal colCount_Y_inout = Y_inout->domain()->dim();

    Teuchos::ArrayRCP<ComplexValueType> complexArray_Y_inout;
    RTOpPack::SubMultiVectorView<ComplexValueType> complexView_Y_inout;
    Teuchos::RCP<Thyra::MultiVectorBase<ComplexValueType> > complex_Y_inout;

    Thyra::DetachedMultiVectorView<ValueType> view_Y_inout(
                Teuchos::rcpFromRef(*Y_inout));
    if (view_Y_inout.leadingDim() == rowCount_Y_inout) {
        // contiguous vector
        ComplexValueType* complexData_Y_inout =
                reinterpret_cast<ComplexValueType*>(view_Y_inout.values());
        complexArray_Y_inout = Teuchos::arcp(complexData_Y_inout,
                                             0, // lowerOffset
                                             (rowCount_Y_inout / 2) * colCount_Y_inout,
                                             false); // doesn't own
        complexView_Y_inout.initialize(0, // globalOffset
                                       rowCount_Y_inout / 2,
                                       0, // colOffset
                                       colCount_Y_inout,
                                       complexArray_Y_inout,
                                       rowCount_Y_inout / 2); // leadingDim
        complex_Y_inout = Thyra::createMembersView(m_complexOperator->range(),
                                                   complexView_Y_inout);
    } else
        throw std::runtime_error("RealWrapperOfComplexThyraLinearOperator::"
                                 "applyImpl(): discontiguous multivectors are not "
                                 "supported yet");

    return m_complexOperator->apply(M_trans, *complex_X_in, complex_Y_inout.ptr(),
                                    alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(
        RealWrapperOfComplexThyraLinearOperator);

} // namespace Bempp

#endif // WITH_TRILINOS
