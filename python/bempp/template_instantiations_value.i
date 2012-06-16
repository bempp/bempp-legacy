// Macros for classes templated on value type

%define BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(CLASS)
template <typename ValueType> class CLASS;
%enddef

%define BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_VALUE(CLASS)

%template(CLASS ## _float32)
    CLASS<float>;
%template(CLASS ## _complex64)
    CLASS<std::complex<float> >;

%template(CLASS ## _float64)
    CLASS<double>;
%template(CLASS ## _complex128)
    CLASS<std::complex<double> >;

%enddef // BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_VALUE

%define BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_VALUE(CLASS)

%extend CLASS<float>
{
    %pythonappend CLASS
    %{ 
        self._valueType = "float32"
    %}
}

%extend CLASS<double>
{
    %pythonappend CLASS
    %{ 
        self._valueType = "float64"
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythonappend CLASS
    %{ 
        self._valueType = "complex64"
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythonappend CLASS
    %{ 
        self._valueType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_VALUE

%define BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_VALUE(CLASS)
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(CLASS);
BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_VALUE(CLASS);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_VALUE(CLASS);
%enddef // BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_VALUE

%pythoncode 
%{

def constructObjectTemplatedOnValue(className, valueType,
                                    *args, **kwargs):
    fullName = className + "_" + valueType
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
