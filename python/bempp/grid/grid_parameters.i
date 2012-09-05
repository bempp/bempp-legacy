%{
#include "grid/grid_parameters.hpp"
%}

namespace Bempp
{

%typemap(out) (GridParameters::Topology)
{
  Bempp::GridParameters::Topology topology = $1;
  if (topology == Bempp::GridParameters::LINEAR)
    $result = PyString_FromString("linear");
  else if (topology == Bempp::GridParameters::TRIANGULAR)
    $result = PyString_FromString("triangular");
  else if (topology == Bempp::GridParameters::QUADRILATERAL)
    $result = PyString_FromString("quadrilateral");
  else {
    PyErr_SetString(PyExc_TypeError, "in method $symname', topology type not supported");
    SWIG_fail;
  }
}

%extend GridParameters
{
%feature("autodoc", "topology -> string") topology;
}

} // namespace Bempp

%include "grid/grid_parameters.hpp"
