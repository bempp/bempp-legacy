#include "context.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "default_local_assembler_factory_for_operators_on_surfaces.hpp"

#include <boost/make_shared.hpp>

namespace Bempp
{

//template <typename BasisFunctionType, typename ResultType>
//shared_ptr<const Context<BasisFunctionType, ResultType> >
//Context<BasisFunctionType, ResultType>::m_defaultContext(
//        boost::make_shared<Context<BasisFunctionType, ResultType>(
//            AssemblyOptions(),
//            boost::make_shared<DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<
//            BasisFunctionType, ResultType> >()));

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
        const AssemblyOptions& assemblyOptions,
        const shared_ptr<LocalAssemblerFactory>& localAssemblerFactory) :
    m_assemblyOptions(assemblyOptions),
    m_localAssemblerFactory(localAssemblerFactory)
{
    if (localAssemblerFactory.get() == 0)
        throw std::invalid_argument("Context::Context(): "
                                    "localAssemblerFactory cannot be null");
}

//template <typename BasisFunctionType, typename ResultType>
//void Context<BasisFunctionType, ResultType>::setDefaultContext(
//        const shared_ptr<const Context>& newContext)
//{
//    if (newContext.get() == 0)
//        throw std::invalid_argument("Context::setDefaultContext(): "
//                                    "new default context cannot be null");
//    m_defaultContext = newContext;
//}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp
