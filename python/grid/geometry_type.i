%{
#include "grid/geometry_type.hpp"
%}

// Handle the enum Dune::GeometryType::BasicType like a string
%typemap(typecheck) Dune::GeometryType::BasicType
{
    $1 = PyString_Check($input);
}
%typemap(in) Dune::GeometryType::BasicType 
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument "
            "$argnum: expected a string"); 
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "simplex")
        $1 = Dune::GeometryType::simplex;
    else if (s == "cube")
        $1 = Dune::GeometryType::cube;
    else if (s == "pyramid")
        $1 = Dune::GeometryType::pyramid;
    else if (s == "prism")
        $1 = Dune::GeometryType::prism;
    else if (s == "extended")
        $1 = Dune::GeometryType::extended;
    else if (s == "none")
        $1 = Dune::GeometryType::none;
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', argument "
            "$argnum: expected one of 'simplex', 'cube', 'pyramid', 'prism', "
            "'extended', or 'none'");        
        SWIG_fail;
    }
}

namespace Dune {

%extend GeometryType
{
    GeometryType(unsigned int d) {
        if (d < 0 || d > 1) {
            throw std::invalid_argument("this constructor works only for dimensions 0 and 1");        
        }
        return new Dune::GeometryType(d);
    }
    %ignore GeometryType(unsigned int, unsigned int);
    %ignore GeometryType(unsigned int);
    %ignore GeometryType(int);

    %ignore basicType;
    %ignore BasicType;
    %ignore Binary;
}

%ignore operator<<;
    
} // namespace Dune

#define DUNE_DEPRECATED
%include <dune/common/geometrytype.hh>
