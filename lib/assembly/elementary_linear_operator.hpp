#ifndef bempp_elementary_linear_operator_hpp
#define bempp_elementary_linear_operator_hpp

#include "linear_operator.hpp"

#include <vector>

namespace Fiber
{

template <typename ValueType> class LocalAssemblerForOperators;
template <typename ValueType> class RawGridGeometry;
template <typename ValueType> class Basis;
template <typename CoordinateType, typename IndexType> class OpenClHandler;

} // namespace Fiber

namespace Bempp
{

template <typename ValueType>
class ElementaryLinearOperator : public LinearOperator<ValueType>
{
public:
    typedef typename LinearOperator<ValueType>::LocalAssemblerFactory
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForOperators<ValueType> LocalAssembler;

    explicit ElementaryLinearOperator(ValueType multiplier = 1.) :
        m_multiplier(multiplier) {
    }

    /** \brief Using a specified factory, construct a local assembler suitable
       for this operator. */
    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const GeometryFactory& geometryFactory,
            const Fiber::RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Fiber::Basis<ValueType>*>& testBases,
            const std::vector<const Fiber::Basis<ValueType>*>& trialBases,
            const Fiber::OpenClHandler<ValueType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** \brief Assemble the operator's weak form using a specified local assembler.

      This function is not intended to be called directly by the user. */
    virtual std::auto_ptr<DiscreteLinearOperator<ValueType> >
    assembleWeakFormInternal(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const = 0;

    /** \brief Multiply the operator in-place by a scalar.

      This method affects the results of subsequent calls to assembleWeakForm()
      and assembleOperator(). */
    void scale(ValueType multiplier) {
        m_multiplier = multiplier;
    }

    /** \brief Return the current value of the scalar by which this operator is multiplied. */
    ValueType multiplier() const {
        return m_multiplier;
    }

private:
    ValueType m_multiplier;
};

}

#endif // bempp_elementary_linear_operator_hpp
