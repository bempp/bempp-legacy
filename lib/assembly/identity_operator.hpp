#ifndef bempp_identity_operator_hpp
#define bempp_identity_operator_hpp

#include "linear_operator.hpp"
#include "../fiber/scalar_function_value.hpp"

namespace Fiber
{

template <typename ValueType> class LocalAssemblerForIdentityOperator;

} // namespace Fiber

namespace Bempp
{

template <typename ValueType>
class IdentityOperator : public LinearOperator<ValueType>
{
public:
    typedef typename LinearOperator<ValueType>::LocalAssemblerFactory
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForIdentityOperator<ValueType> LocalAssembler;

    virtual int trialComponentCount() const { return 1; }

    virtual int testComponentCount() const { return 1; }

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
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInDenseMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInSparseMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    Fiber::ScalarFunctionValue<ValueType> m_expression;
};

} // namespace Bempp

#endif
