%{
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/grid_parameters.hpp"
%}

%include "grid_docstrings.i"

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

    
%extend Grid 
{
    %pythonappend leafView %{
        val._parentGrid = self
    %}

    %pythonappend levelView %{
        val._parentGrid = self
    %}

    %pythonappend globalIdSet %{
        val._parentGrid = self
    %}

    // this function is only for internal use
    %ignore elementGeometryFactory;

}

} // namespace Bempp

%include "grid/grid.hpp"

%clear GridParameters::Topology;
