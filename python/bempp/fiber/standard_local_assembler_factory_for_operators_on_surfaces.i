namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces :
    public LocalAssemblerFactory<
        BasisFunctionType, ResultType, GeometryFactory>
{
};

} // namespace Fiber

namespace Fiber
{
BEMPP_PYTHON_INSTANTIATE_ANONYMOUSLY_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
StandardLocalAssemblerFactoryForOperatorsOnSurfaces);
}
