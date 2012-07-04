#include "config_trilinos.hpp"

#ifndef bempp_real_wrapper_of_complex_thyra_preconditioner_hpp
#define bempp_real_wrapper_of_complex_thyra_preconditioner_hpp

#include "../common/common.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_PreconditionerBase.hpp>

namespace Bempp
{

template <typename ValueType>
class RealWrapperOfComplexThyraPreconditioner
        : public Thyra::PreconditionerBase<ValueType>
{
public:
    typedef std::complex<ValueType> ComplexValueType;
    typedef Thyra::PreconditionerBase<std::complex<ValueType> > ComplexPreconditioner;

    RealWrapperOfComplexThyraPreconditioner(
            const Teuchos::RCP<const ComplexPreconditioner>& complexPreconditioner);

    /** \brief Return if the underlying left preconditioner operator is
     * const-only or allows non-const access.
     */
    virtual bool isLeftPrecOpConst() const;

    /** \brief Return a non-const left preconditioner linear operator if one is
     * designed or targeted to be applied on the left.
     *
     * <b>Preconditions:</b><ul>
     * <li>[<tt>isLeftPrecOpConst()==true</tt>] <tt>getLeftPrecOp().get()==NULL</tt>
     * </ul>
     */
    virtual Teuchos::RCP<Thyra::LinearOpBase<ValueType> > getNonconstLeftPrecOp();

    /** \brief Return a const left preconditioner linear operator if one is
     * designed or targeted to be applied on the left.
     */
    virtual Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > getLeftPrecOp() const;

    /** \brief Return if the underlying right preconditioner operator is
     * const-only or allows non-const access.
     */
    virtual bool isRightPrecOpConst() const;

    /** \brief Return a non-const right preconditioner linear operator if one is
     * designed or targeted to be applied on the right.
     *
     * <b>Preconditions:</b><ul>
     * <li>[<tt>isRightPrecOpConst()==true</tt>] <tt>getRightPrecOp().get()==NULL</tt>
     * </ul>
     */
    virtual Teuchos::RCP<Thyra::LinearOpBase<ValueType> > getNonconstRightPrecOp();

    /** \brief Return a const right preconditioner linear operator if one is
     * designed or targeted to be applied on the right.
     */
    virtual Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > getRightPrecOp() const;

    /** \brief Return if the underlying unspecified preconditioner operator is
     * const-only or allows non-const access.
     */
    virtual bool isUnspecifiedPrecOpConst() const;

    /** \brief Return a non-const generic preconditioner linear operator that is
     * not designed or targeted to be applied on the left or on the right.
     */
    virtual Teuchos::RCP<Thyra::LinearOpBase<ValueType> > getNonconstUnspecifiedPrecOp();

    /** \brief Return a const generic preconditioner linear operator that is not
     * designed or targeted to be applied on the left or on the right.
     *
     * <b>Preconditions:</b><ul>
     * <li>[<tt>isUnspecifiedPrecOpConst()==true</tt>] <tt>getUnspecifiedPrecOp().get()==NULL</tt>
     * </ul>
     */
    virtual Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > getUnspecifiedPrecOp() const;

private:
    Teuchos::RCP<const ComplexPreconditioner> m_complexPreconditioner;
};


} // namespace Bempp

#endif // WITH_TRILINOS

#endif
