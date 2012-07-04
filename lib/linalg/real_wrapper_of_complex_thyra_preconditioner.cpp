#include "real_wrapper_of_complex_thyra_preconditioner.hpp"

#ifdef WITH_TRILINOS

#include "real_wrapper_of_complex_thyra_linear_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename ValueType>
RealWrapperOfComplexThyraPreconditioner<ValueType>::RealWrapperOfComplexThyraPreconditioner(
        const Teuchos::RCP<const ComplexPreconditioner>& complexPreconditioner) :
    m_complexPreconditioner(complexPreconditioner)
{
    if (m_complexPreconditioner.is_null())
        throw std::invalid_argument(
                "RealWrapperOfComplexThyraPreconditioner::"
                "RealWrapperOfComplexThyraPreconditioner(): "
                "argument must not be null");
}

template <typename ValueType>
bool
RealWrapperOfComplexThyraPreconditioner<ValueType>::isLeftPrecOpConst() const
{
    return false;
}

template <typename ValueType>
Teuchos::RCP<Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getNonconstLeftPrecOp()
{
    throw std::runtime_error(
                "RealWrapperOfComplexThyraPreconditioner::getNonconstLeftPrecOp(): "
                "preconditioner is constant");
}

template <typename ValueType>
Teuchos::RCP<const Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getLeftPrecOp() const
{
    Teuchos::RCP<const Thyra::LinearOpBase<ComplexValueType> > complexPrecOp =
            m_complexPreconditioner->getLeftPrecOp();
    if (complexPrecOp.is_null())
        return Teuchos::null;
    else
        return Teuchos::rcp<const Thyra::LinearOpBase<ValueType> >(
                    new RealWrapperOfComplexThyraLinearOperator<ValueType>(
                        complexPrecOp));
}

template <typename ValueType>
bool
RealWrapperOfComplexThyraPreconditioner<ValueType>::isRightPrecOpConst() const
{
    return false;
}

template <typename ValueType>
Teuchos::RCP<Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getNonconstRightPrecOp()
{
    throw std::runtime_error(
                "RealWrapperOfComplexThyraPreconditioner::getNonconstRightPrecOp(): "
                "preconditioner is constant");
}

template <typename ValueType>
Teuchos::RCP<const Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getRightPrecOp() const
{
    Teuchos::RCP<const Thyra::LinearOpBase<ComplexValueType> > complexPrecOp =
            m_complexPreconditioner->getRightPrecOp();
    if (complexPrecOp.is_null())
        return Teuchos::null;
    else
        return Teuchos::rcp<const Thyra::LinearOpBase<ValueType> >(
                    new RealWrapperOfComplexThyraLinearOperator<ValueType>(
                        complexPrecOp));

}

template <typename ValueType>
bool
RealWrapperOfComplexThyraPreconditioner<ValueType>::isUnspecifiedPrecOpConst() const
{
    return false;
}

template <typename ValueType>
Teuchos::RCP<Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getNonconstUnspecifiedPrecOp()
{
    throw std::runtime_error(
                "RealWrapperOfComplexThyraPreconditioner::getNonconstUnspecifiedPrecOp(): "
                "preconditioner is constant");
}

template <typename ValueType>
Teuchos::RCP<const Thyra::LinearOpBase<ValueType> >
RealWrapperOfComplexThyraPreconditioner<ValueType>::getUnspecifiedPrecOp() const
{
    Teuchos::RCP<const Thyra::LinearOpBase<ComplexValueType> > complexPrecOp =
            m_complexPreconditioner->getUnspecifiedPrecOp();
    if (complexPrecOp.is_null())
        return Teuchos::null;
    else
        return Teuchos::rcp<const Thyra::LinearOpBase<ValueType> >(
                    new RealWrapperOfComplexThyraLinearOperator<ValueType>(
                        complexPrecOp));

}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(
        RealWrapperOfComplexThyraPreconditioner);

} // namspace Bempp

#endif // WITH_TRILINOS
