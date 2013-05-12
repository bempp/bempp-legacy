%{
#include "assembly/aca_options.hpp"
%}

namespace Bempp
{

%feature("autodoc", "eps -> float") AcaOptions::eps;
%feature("autodoc", "eta -> float") AcaOptions::eta;
%feature("autodoc", "mode -> \"global_assembly\", \"local_assembly\" or \"hybrid_assembly\"")
     AcaOptions::mode;
%feature("autodoc", "reactionToUnsupportedMode -> \"ignore\", \"warning\" or \"error\"")
     AcaOptions::reactionToUnsupportedMode;
%feature("autodoc", "maximumBlockSize -> int") AcaOptions::maximumBlockSize;
%feature("autodoc", "maximumRank -> int") AcaOptions::maximumRank;
%feature("autodoc", "minimumBlockSize -> int") AcaOptions::minimumBlockSize;
%feature("autodoc", "outputFname -> string") AcaOptions::outputFname;
%feature("autodoc", "outputPostscript -> bool") AcaOptions::outputPostscript;
%feature("autodoc", "recompress -> bool") AcaOptions::recompress;
%feature("autodoc", "scaling -> float") AcaOptions::scaling;
%feature("autodoc", "firstClusterIndex -> int") AcaOptions::firstClusterIndex;
%feature("autodoc", "globalAssemblyBeforeCompression (deprecated) -> bool")
     AcaOptions::globalAssemblyBeforeCompression;

// Handle the enum AcaOptions::AcaAssemblyMode like a string
%typemap(typecheck) AcaOptions::AcaAssemblyMode
{
    $1 = PyString_Check($input);
}
%typemap(in) AcaOptions::AcaAssemblyMode
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, 
                        "in method '$symname', argument $argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "global_assembly")
        $1 = Bempp::AcaOptions::GLOBAL_ASSEMBLY;
    else if (s == "local_assembly")
        $1 = Bempp::AcaOptions::LOCAL_ASSEMBLY;
    else if (s == "hybrid_assembly")
        $1 = Bempp::AcaOptions::HYBRID_ASSEMBLY;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of 'global_assembly', "
                        "'local_assembly' or 'hybrid_assembly'");
        SWIG_fail;
    }
}

%typemap(out) AcaOptions::AcaAssemblyMode
{
    if ($1 == Bempp::AcaOptions::GLOBAL_ASSEMBLY)
        $result = PyString_FromString("global_assembly");
    else if ($1 == Bempp::AcaOptions::LOCAL_ASSEMBLY)
        $result = PyString_FromString("local_assembly");
    else if ($1 == Bempp::AcaOptions::HYBRID_ASSEMBLY)
        $result = PyString_FromString("hybrid_assembly");
    else
    {
        PyErr_SetString(PyExc_ValueError, 
                        "in method '$symname', unknown ACA assembly mode.");
        SWIG_fail;
    }
}

// Handle the enum AcaOptions::ReactionToUnsupportedMode like a string
%typemap(typecheck) AcaOptions::ReactionToUnsupportedMode
{
    $1 = PyString_Check($input);
}
%typemap(in) AcaOptions::ReactionToUnsupportedMode
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, 
                        "in method '$symname', argument $argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "ignore")
        $1 = Bempp::AcaOptions::IGNORE;
    else if (s == "warning")
        $1 = Bempp::AcaOptions::WARNING;
    else if (s == "error")
        $1 = Bempp::AcaOptions::ERROR;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of 'ignore', "
                        "'warning' or 'error'");
        SWIG_fail;
    }
}

%typemap(out) AcaOptions::ReactionToUnsupportedMode
{
    if ($1 == Bempp::AcaOptions::IGNORE)
        $result = PyString_FromString("ignore");
    else if ($1 == Bempp::AcaOptions::WARNING)
        $result = PyString_FromString("warning");
    else if ($1 == Bempp::AcaOptions::ERROR)
        $result = PyString_FromString("error");
    else
    {
        PyErr_SetString(PyExc_ValueError, 
                        "in method '$symname', unknown reaction to unsupported mode.");
        SWIG_fail;
    }
}

} // namespace Bempp

%include "assembly/aca_options.hpp"
