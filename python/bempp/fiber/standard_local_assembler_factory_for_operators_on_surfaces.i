namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces :
    public LocalAssemblerFactory<
        BasisFunctionType, ResultType, GeometryFactory>
{
//public:
//    const AccuracyOptions& accuracyOptions();
};

} // namespace Fiber

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
StandardLocalAssemblerFactoryForOperatorsOnSurfaces);
}
