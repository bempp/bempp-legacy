#ifndef fiber_integration_manager_factory_hpp
#define fiber_integration_manager_factory_hpp

#include <armadillo>
#include <memory>

namespace Fiber
{

template <typename ValueType> class Expression;
template <typename ValueType> class Kernel;
template <typename ValueType> class IntegrationManager;

template <typename ValueType, typename GeometryFactory>
class IntegrationManagerFactory
{
public:
    virtual ~IntegrationManagerFactory() {}

    virtual std::auto_ptr<IntegrationManager<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const arma::Mat<ValueType>& vertices,
            const arma::Mat<int>& elementCornerIndices,
            const arma::Mat<char>& auxData,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression) const = 0;
};

} // namespace Fiber

#endif
