namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class LocalAssemblerFactory
{
public:
    virtual ~LocalAssemblerFactory() = 0;
};

} // namespace Fiber

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
LocalAssemblerFactory);
}
