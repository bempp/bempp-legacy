#ifdef WITH_TRILINOS

%{
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <sstream>
%}

// Temporary solution: convert parameter lists to XML strings.
// In future we may want to convert them to Python dictionaries.

%typemap(typecheck) const Teuchos::RCP<Teuchos::ParameterList>&
{
    $1 = PyString_Check($input);
}
%typemap(in) const Teuchos::RCP<Teuchos::ParameterList>&
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument "
            "$argnum: expected a string");
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    $1 = new Teuchos::RCP<Teuchos::ParameterList>(Teuchos::getParametersFromXmlString(s));
}
%typemap(freearg) const Teuchos::RCP<Teuchos::ParameterList>&
{
    delete $1;
}

%typemap(out) Teuchos::RCP<Teuchos::ParameterList>
{
    std::ostringstream out;
    Teuchos::writeParameterListToXmlOStream(*$1, out);
    $result = PyString_FromString(out.str().c_str());
}

#endif
