%{
#include "grid/vtk_writer.hpp"
%}

namespace Dune 
{
namespace VTK 
{

/** \bug Improve error reporting */
%typemap(in) OutputType 
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "Expected a string");        
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "ascii")
        $1 = ascii;
    else if (s == "base64")
        $1 = base64;
    else if (s == "appendedraw")
        $1 = appendedraw;
    else if (s == "appendedbase64")
        $1 = appendedbase64;
    else
        SWIG_fail;
}

} // namespace VTK
} // namespace Dune
    
namespace Bempp 
{
    
%apply const arma::Mat<double>& IN_MAT_OUT_WRAPPERS {
    const arma::Mat<double>& data
}; 

%extend VtkWriter 
{
    // Store references to the Numpy and Armadillo array objects 
    // holding data to be exported to the VTK file
    %pythonappend addCellData %{
        try:
            self._python_data.append(val[0])
        except AttributeError:
            self._python_data = [val[0]]
        try:
            self._c_data.append(val[1])
        except AttributeError:
            self._c_data = [val[1]]
        val = None
    %}

    // Store references to the Numpy and Armadillo array objects 
    // holding data to be exported to the VTK file
    %pythonappend addVertexData %{
        try:
            self._python_data.append(val[0])
        except AttributeError:
            self._python_data = [val[0]]
        try:
            self._c_data.append(val[1])
        except AttributeError:
            self._c_data = [val[1]]
        val = None
    %}

    // No references to Numpy or Armadillo array objects 
    // need to be stored any longer
    %pythonappend clear %{
        self._python_data = None
        self._c_data = None
    %}

    // Reference to the parent grid view, stored to prevent it from
    // being garbage-collected while this VTK writer is alive
    %pythoncode %{
        def parentGridView(self): return self._parentGridView
        parentGridView = property(parentGridView, None, None, "Parent grid view")
    %}
}

} // namespace Bempp

%include "grid/vtk_writer.hpp"

namespace Bempp {
%clear const arma::Mat<double>& data;
}
