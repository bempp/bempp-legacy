%{
#include "assembly/assembly_options.hpp"
%}

namespace Bempp
{

// Handle the enum AssemblyOptions::Value like a string
%typemap(typecheck) AssemblyOptions::Value
{
    $1 = PyString_Check($input);
}
%typemap(in) AssemblyOptions::Value
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,
                        "in method '$symname', argument $argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "auto")
        $1 = Bempp::AssemblyOptions::AUTO;
    else if (s == "yes")
        $1 = Bempp::AssemblyOptions::YES;
    else if (s == "no")
        $1 = Bempp::AssemblyOptions::NO;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of 'auto', "
                        "'yes' or 'no'");
        SWIG_fail;
    }
}

%typemap(out) AssemblyOptions::Value
{
    if ($1 == Bempp::AssemblyOptions::AUTO)
        $result = PyString_FromString("auto");
    else if ($1 == Bempp::AssemblyOptions::YES)
        $result = PyString_FromString("yes");
    else if ($1 == Bempp::AssemblyOptions::NO)
        $result = PyString_FromString("no");
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', invalid value.");
        SWIG_fail;
    }
}

%extend AssemblyOptions
{
    %ignore switchToTbb;
    %feature("compactdefaultargs") enableSingularIntegralCaching;
    %feature("compactdefaultargs") enableSparseStorageOfMassMatrices;
    %feature("compactdefaultargs") enableJointAssembly;
    %feature("compactdefaultargs") enableBlasInQuadrature;
}

} // namespace Bempp

%include "assembly/assembly_options.hpp"
