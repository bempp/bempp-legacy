// Macros for classes templated on result type

%define BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_RESULT(CLASS)
template <typename ResultType> class CLASS;
%enddef

%define BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(CLASS)

%template(CLASS ## _float32)
    CLASS<float>;
%template(CLASS ## _complex64)
    CLASS<std::complex<float> >;

%template(CLASS ## _float64)
    CLASS<double>;
%template(CLASS ## _complex128)
    CLASS<std::complex<double> >;

%enddef // BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT

%define BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_RESULT(CLASS)

%extend CLASS<float>
{
    %pythonappend CLASS
    %{ 
        self._resultType = "float32"
    %}
}

%extend CLASS<double>
{
    %pythonappend CLASS
    %{ 
        self._resultType = "float64"
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythonappend CLASS
    %{ 
        self._resultType = "complex64"
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythonappend CLASS
    %{ 
        self._resultType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_RESULT

%define BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_RESULT(CLASS)
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_RESULT(CLASS);
BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_RESULT(CLASS);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(CLASS);
%enddef // BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_RESULT

%pythoncode 
%{

def constructObjectTemplatedOnResult(className, resultType,
                                    *args, **kwargs):
    fullName = className + "_" + resultType
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
