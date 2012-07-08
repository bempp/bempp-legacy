namespace Fiber
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(LocalAssemblerFactory);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(LocalAssemblerFactory);
}

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class LocalAssemblerFactory
{
/* public: */
/*     virtual ~LocalAssemblerFactory() = 0; */
};

} // namespace Fiber

%shared_ptr(Fiber::LocalAssemblerFactory<
            float, float, Bempp::GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<
            float, std::complex<float>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<
            std::complex<float>, std::complex<float>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<
            double, double, Bempp::GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<
            double, std::complex<double>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::LocalAssemblerFactory<
            std::complex<double>, std::complex<double>, Bempp::GeometryFactory>);

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
LocalAssemblerFactory);
}
