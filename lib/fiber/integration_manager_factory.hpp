#ifndef fiber_integration_manager_factory_hpp
#define fiber_integration_manager_factory_hpp

#include <memory>

namespace Fiber
{

template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;
template <typename ValueType, typename GeometryImp> class IntegrationManager;

template <typename ValueType, typename GeometryImp>
class IntegrationManagerFactory
{
public:
    virtual ~IntegrationManagerFactory() {}

    virtual std::auto_ptr<IntegrationManager<ValueType, GeometryImp> > make(
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression) const = 0;
};

} // namespace Fiber

#endif
