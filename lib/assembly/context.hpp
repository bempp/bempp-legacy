#ifndef bempp_context_hpp
#define bempp_context_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "assembly_options.hpp"
#include "discrete_boundary_operator_cache.hpp"

namespace Bempp
{

class GeometryFactory;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class AbstractBoundaryOperator;

template <typename BasisFunctionType, typename ResultType>
class Context
{
public:
    /** \brief Type of the appropriate instantiation of Fiber::LocalAssemblerFactory. */
    typedef Fiber::LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory>
    LocalAssemblerFactory;

    Context(const AssemblyOptions& assemblyOptions,
            const shared_ptr<LocalAssemblerFactory>& localAssemblerFactory);

//    static shared_ptr<const Context>& defaultContext() {
//        return m_defaultContext;
//    }

//    static void setDefaultContext(const shared_ptr<const Context>& newContext);

    shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    getWeakForm(const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op) const {
        return m_cache.getWeakForm(*this, op);
    }

private:
//    static shared_ptr<const Context> m_defaultContext;

    AssemblyOptions m_assemblyOptions;
    shared_ptr<LocalAssemblerFactory> m_localAssemblerFactory;

    mutable DiscreteBoundaryOperatorCache<BasisFunctionType, ResultType> m_cache;
};

} // namespace Bempp

#endif
