#ifndef fiber_standard_integration_manager_factory_2d_hpp
#define fiber_standard_integration_manager_factory_2d_hpp

#include "integration_manager_factory.hpp"
#include "standard_integration_manager_2d.hpp"
#include "opencl_options.hpp"

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
class StandardIntegrationManagerFactory2D :
        public IntegrationManagerFactory<ValueType, GeometryFactory>
{
public:
    StandardIntegrationManagerFactory2D(const OpenClOptions& openClOptions) :
        m_openClOptions(openClOptions)
    {}

    virtual std::auto_ptr<IntegrationManager<ValueType> > make(
            const GeometryFactory& geometryFactory,
            const arma::Mat<ValueType>& vertices,
            const arma::Mat<int>& elementCornerIndices,
            const arma::Mat<char>& auxData,
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression) const
    {
        return std::auto_ptr<IntegrationManager<ValueType> >(
                    new StandardIntegrationManager2D<ValueType, GeometryFactory>(
                        geometryFactory, vertices, elementCornerIndices, auxData,
                        testExpression, kernel, trialExpression,
                        m_openClOptions));
    }

private:
    OpenClOptions m_openClOptions;
};

} // namespace Fiber

#endif
