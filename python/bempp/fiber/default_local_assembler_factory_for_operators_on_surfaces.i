namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class DefaultLocalAssemblerFactoryForOperatorsOnSurfaces :
    public LocalAssemblerFactory<
        BasisFunctionType, ResultType, GeometryFactory>
{
//public:
//    const AccuracyOptions& accuracyOptions();
};

} // namespace Fiber

%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<float, float, Bempp::GeometryFactory>);
%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<float, std::complex<float>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<float>, std::complex<float>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<double, double, Bempp::GeometryFactory>);
%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<double, std::complex<double>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<std::complex<double>, std::complex<double>, Bempp::GeometryFactory >);

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
DefaultLocalAssemblerFactoryForOperatorsOnSurfaces);
}
