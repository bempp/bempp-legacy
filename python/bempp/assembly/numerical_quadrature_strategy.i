%{
#include "assembly/numerical_quadrature_strategy.hpp"
%}

%define BEMPP_NUMERICAL_QUADRATURE_STRATEGY_FACTORY(BASIS,RESULT,PY_BASIS,PY_RESULT)
%inline %{

namespace Bempp
{

boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >
numericalQuadratureStrategy_## PY_BASIS ## _ ## PY_RESULT(
    const Fiber::AccuracyOptionsEx& accuracyOptions)
{
    return boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >(
       new NumericalQuadratureStrategy< BASIS , RESULT >(accuracyOptions));
}

boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >
numericalQuadratureStrategy_## PY_BASIS ## _ ## PY_RESULT(
    const Fiber::AccuracyOptions& accuracyOptions)
{
    return boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >(
       new NumericalQuadratureStrategy< BASIS , RESULT >(accuracyOptions));
}

} // namespace Bempp

%}
%enddef

BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(BEMPP_NUMERICAL_QUADRATURE_STRATEGY_FACTORY)

