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

template <typename ValueType> class LocalAssemblerForIntegralOperators;
template <typename ValueType> class LocalAssemblerForIdentityOperator;
template <typename ValueType> class LocalAssemblerForSourceTerms;

template <typename ValueType, typename GeometryFactory>
class LocalAssemblerFactory
{
public:
    virtual ~LocalAssemblerFactory() {}

    /** @name Local assemblers for integral operators
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for an integral operator. */
    virtual std::auto_ptr<LocalAssemblerForIntegralOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an integral operator.

        Used also for evaluation of the identity operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForIntegralOperators<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** @}
        @name Local assemblers for the identity operator
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for the identity operator. */
    virtual std::auto_ptr<LocalAssemblerForIdentityOperator<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an identity operator.

        Used also for evaluation of the identity operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForIdentityOperator<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& trialExpression,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** @}
        @name Local assemblers for source terms
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for a source term. */
    virtual std::auto_ptr<LocalAssemblerForSourceTerms<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const Expression<ValueType>& testExpression,
            const Function<ValueType>& function,
            const OpenClHandler<ValueType,int>& openClHandler) const = 0;

    /** @} */
};

} // namespace Fiber

#endif
