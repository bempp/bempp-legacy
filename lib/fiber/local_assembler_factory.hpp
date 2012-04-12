#ifndef fiber_integration_manager_factory_hpp
#define fiber_integration_manager_factory_hpp

#include <armadillo>
#include <memory>

namespace Fiber
{
template <typename ValueType, typename IndexType> class OpenClHandler;

template <typename ValueType> class Basis;
template <typename ValueType> class Expression;
template <typename ValueType> class Function;
template <typename ValueType> class Kernel;
template <typename ValueType> class RawGridGeometry;

template <typename ValueType> class LocalAssemblerForOperators;
template <typename ValueType> class LocalAssemblerForGridFunctions;
template <typename ValueType> class EvaluatorForIntegralOperators;

template <typename ValueType, typename GeometryFactory>
class LocalAssemblerFactory
{
public:
    virtual ~LocalAssemblerFactory() {}

    /** @name Local assemblers for integral operators
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for an integral operator. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            ValueType multiplier,
            const OpenClHandler<ValueType,int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an integral operator.

        Used also for evaluation of the identity operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            ValueType multiplier,
            const OpenClHandler<ValueType,int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** @}
        @name Local assemblers for the identity operator
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for the identity operator. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Expression<ValueType>& trialExpression,
            ValueType multiplier,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an identity operator.

        Used also for evaluation of the identity operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& trialExpression,
            ValueType multiplier,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** @}
        @name Local assemblers for grid functions
        @{ */

    /** \brief Allocate a local assembler for calculations of the projections
      of functions from a given space on a Fiber::Function. */
    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const Expression<ValueType>& testExpression,
            const Function<ValueType>& function,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** @}
        @name Evaluators for integral operators
        @{ */

    /** \brief Allocate an evaluator for an integral operator applied to a
      grid function. */
    virtual std::auto_ptr<EvaluatorForIntegralOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const std::vector<std::vector<ValueType> >& argumentLocalCoefficients,
            ValueType multiplier,
            const OpenClHandler<ValueType, int>& openClHandler) const = 0;

    /** @} */
};

} // namespace Fiber

#endif
