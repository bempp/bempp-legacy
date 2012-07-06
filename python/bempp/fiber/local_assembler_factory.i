%shared_ptr(Fiber::LocalAssemblerFactory<float, float, GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<float, std::complex<float> , GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<std::complex<float>, std::complex<float> , GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<double, double, GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<double, std::complex<double> , GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<std::complex<double>, std::complex<double> , GeometryFactory>);

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
