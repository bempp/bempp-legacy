%{
#include "grid/geometry_factory.hpp"
#include "fiber/quadrature_strategy.hpp"
%}

namespace Fiber
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategyBase);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategyBase);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategy);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(QuadratureStrategy);

 %warnfilter(520) QuadratureStrategy;

%extend QuadratureStrategyBase {
%ignore makeAssemblerForIntegralOperators;
%ignore makeAssemblerForIdentityOperators;
%ignore makeAssemblerForLocalOperators;
%ignore makeAssemblerForGridFunctions;
%ignore makeEvaluatorForIntegralOperators;
}

}


#define shared_ptr boost::shared_ptr
%include "fiber/quadrature_strategy.hpp"
#undef shared_ptr

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
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
QuadratureStrategyBase);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(
QuadratureStrategy);
}
