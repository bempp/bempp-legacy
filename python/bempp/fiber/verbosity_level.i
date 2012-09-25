%{
#include "fiber/verbosity_level.hpp"
%}

// Handle the enum Fiber::VerbosityLevel::Level like a string
%typemap(typecheck) Fiber::VerbosityLevel::Level
{
    $1 = PyString_Check($input);
}
%typemap(in) Fiber::VerbosityLevel::Level
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument $argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "default")
        $1 = Fiber::VerbosityLevel::DEFAULT;
    else if (s == "low")
        $1 = Fiber::VerbosityLevel::LOW;
    else if (s == "high")
        $1 = Fiber::VerbosityLevel::HIGH;
    else
    {
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected one of 'default', 'low' or 'high'");
        SWIG_fail;
    }
}

%typemap(out) Fiber::VerbosityLevel::Level
{
    if ($1 == Fiber::VerbosityLevel::DEFAULT)
        $result = PyString_FromString("default");
    else if ($1 == Fiber::VerbosityLevel::LOW)
        $result = PyString_FromString("low");
    else if ($1 == Fiber::VerbosityLevel::HIGH)
        $result = PyString_FromString("high");
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', unknown verbosity level.");
        SWIG_fail;
    }
}

%include "fiber/verbosity_level.hpp"

%apply Fiber::VerbosityLevel::Level { Bempp::VerbosityLevel::Level };
