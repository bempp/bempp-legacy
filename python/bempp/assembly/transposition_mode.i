%{
#include "assembly/transposition_mode.hpp"
%}

namespace Bempp
{

// Handle the enum TranspositionMode like a string
%typemap(typecheck) TranspositionMode
{
    $1 = PyString_Check($input);
}
%typemap(in) TranspositionMode
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError,
                        "in method '$symname', argument $argnum: "
                        "expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "no_transpose" || s == "n" || s == "")
        $1 = Bempp::NO_TRANSPOSE;
    else if (s == "transpose" || s == "t")
        $1 = Bempp::TRANSPOSE;
    else if (s == "conjugate" || s == "c")
        $1 = Bempp::CONJUGATE;
    else if (s == "conjugate_transpose" || s == "ct" || s == "h")
        $1 = Bempp::CONJUGATE_TRANSPOSE;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of '', 'n', 'no_transpose', "
                        "'t', 'transpose', "
                        "'c', 'conjugate', "
                        "'h', 'ct' or 'conjugate_transpose'");
        SWIG_fail;
    }
}

} // namespace Bempp

%include "assembly/transposition_mode.hpp"
