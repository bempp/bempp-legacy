namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class NumericalQuadratureStrategy :
    public QuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory>
{
  public:
  const AccuracyOptions& accuracyOptions();
  virtual ~NumericalQuadratureStrategy();
};

} // namespace Fiber

%shared_ptr(Fiber::NumericalQuadratureStrategy<float, float, Bempp::GeometryFactory>);
%shared_ptr(Fiber::NumericalQuadratureStrategy<float, std::complex<float>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::NumericalQuadratureStrategy<std::complex<float>, std::complex<float>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::NumericalQuadratureStrategy<double, double, Bempp::GeometryFactory>);
%shared_ptr(Fiber::NumericalQuadratureStrategy<double, std::complex<double>, Bempp::GeometryFactory >);
%shared_ptr(Fiber::NumericalQuadratureStrategy<std::complex<double>, std::complex<double>, Bempp::GeometryFactory >);

namespace Fiber
{
BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
NumericalQuadratureStrategy);
}
