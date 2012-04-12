#ifndef bempp_aca_global_assembler_hpp
#define bempp_aca_global_assembler_hpp

#include <memory>
#include <vector>

namespace Fiber
{

template <typename ValueType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

class AssemblyOptions;
template <typename ValueType> class DiscreteLinearOperator;
template <typename ValueType> class Space;

template <typename ValueType>
class AcaGlobalAssembler
{
public:
    typedef DiscreteLinearOperator<ValueType> DiscreteLinOp;
    typedef Fiber::LocalAssemblerForOperators<ValueType> LocalAssembler;

    static std::auto_ptr<DiscreteLinOp> assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const std::vector<LocalAssembler*>& localAssemblers,
            const std::vector<const DiscreteLinOp*>& sparseTermsToAdd,
            const AssemblyOptions& options);

    static std::auto_ptr<DiscreteLinOp> assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& localAssembler,
            const AssemblyOptions& options);
};

} // namespace Bempp

#endif
