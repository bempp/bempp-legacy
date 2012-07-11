namespace Fiber
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategy);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategy);
}

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class QuadratureStrategy
{
/* public: */
/*     virtual ~QuadratureStrategy() = 0; */
};

} // namespace Fiber

%shared_ptr(Fiber::QuadratureStrategy<
            float, float, Bempp::GeometryFactory>);
%shared_ptr(Fiber::QuadratureStrategy<
            float, std::complex<float>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::QuadratureStrategy<
            std::complex<float>, std::complex<float>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::QuadratureStrategy<
            double, double, Bempp::GeometryFactory>);
%shared_ptr(Fiber::QuadratureStrategy<
            double, std::complex<double>, Bempp::GeometryFactory>);
%shared_ptr(Fiber::QuadratureStrategy<
            std::complex<double>, std::complex<double>, Bempp::GeometryFactory>);

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
QuadratureStrategy);
}
