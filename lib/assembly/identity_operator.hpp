#ifndef bempp_identity_operator_hpp
#define bempp_identity_operator_hpp

#include "elementary_linear_operator.hpp"
#include "../fiber/scalar_function_value.hpp"

namespace Fiber
{

template <typename ValueType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

template <typename ValueType>
class IdentityOperator : public ElementaryLinearOperator<ValueType>
{
public:
    typedef typename ElementaryLinearOperator<ValueType>::LocalAssemblerFactory
    LocalAssemblerFactory;
    typedef typename ElementaryLinearOperator<ValueType>::LocalAssembler LocalAssembler;

    virtual int trialComponentCount() const { return 1; }

    virtual int testComponentCount() const { return 1; }

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const GeometryFactory& geometryFactory,
            const Fiber::RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Fiber::Basis<ValueType>*>& testBases,
            const std::vector<const Fiber::Basis<ValueType>*>& trialBases,
            const Fiber::OpenClHandler<ValueType, int>& openClHandler,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInternal(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
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
