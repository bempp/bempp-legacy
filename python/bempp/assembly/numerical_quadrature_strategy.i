%{
#include "assembly/numerical_quadrature_strategy.hpp"
%}

%define NUMERICAL_QUADRATURE_STRATEGY_FACTORY(BASIS,RESULT,PY_BASIS,PY_RESULT)
     %inline %{
     namespace Bempp{
boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >
     createNumericalQuadratureStrategy_## PY_BASIS ## _ ## PY_RESULT (const Fiber::AccuracyOptions& accuracyOptions)
{
  return boost::shared_ptr<Fiber::QuadratureStrategy< BASIS , RESULT , GeometryFactory> >(
                new NumericalQuadratureStrategy< BASIS , RESULT >(accuracyOptions));
}

 }
     %}
%enddef

BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(NUMERICAL_QUADRATURE_STRATEGY_FACTORY)

