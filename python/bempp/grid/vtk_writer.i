%{
#include "grid/vtk_writer.hpp"
%}

%include "vtk_writer_docstrings.i"

// Handle the enum Dune::VTK::OutputType like a string
%typemap(in) VtkWriter::OutputType
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument $argnum: expected a string"); 
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "ascii")
        $1 = VtkWriter::ASCII;
    else if (s == "base64")
        $1 = VtkWriter::BASE_64;
    else if (s == "appendedraw")
        $1 = VtkWriter::APPENDED_RAW;
    else if (s == "appendedbase64")
        $1 = VtkWriter::APPENDED_BASE_64;
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', argument $argnum: expected one of 'ascii', 'base64', 'appendedraw' or 'appendedbase64'");        
        SWIG_fail;
    }
}
    
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

    %feature("compactdefaultargs") write;
    %feature("compactdefaultargs") pwrite;
}

} // namespace Bempp

%include "grid/vtk_writer.hpp"

namespace Bempp {
%clear const arma::Mat<double>& data;
}
