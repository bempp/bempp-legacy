// Macros for classes templated on basis and result

%define BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
template <typename BasisFunctionType, typename ResultType> class CLASS;
%enddef

%define BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)

%template(CLASS ## _float32_float32)
    CLASS<float, float>;
%template(CLASS ## _float32_complex64)
    CLASS<float, std::complex<float> >;
%template(CLASS ## _complex64_complex64)
    CLASS<std::complex<float>, std::complex<float> >;

%template(CLASS ## _float64_float64)
    CLASS<double, double>;
%template(CLASS ## _float64_complex128)
    CLASS<double, std::complex<double> >;
%template(CLASS ## _complex128_complex128)
    CLASS<std::complex<double>, std::complex<double> >

%enddef // BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT

%define BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)

%extend CLASS<float, float>
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "float32"
        self._resultType = "float32"
    %}
}

%extend CLASS<float, std::complex<float> >
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "float32"
        self._resultType = "complex64"
    %}
}

%extend CLASS<std::complex<float>, std::complex<float> >
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "complex64"
        self._resultType = "complex64"
    %}
}

%extend CLASS<double, double>
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "float64"
        self._resultType = "float64"
    %}
}

%extend CLASS<double, std::complex<double> >
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "float64"
        self._resultType = "complex128"
    %}
}

%extend CLASS<std::complex<double>, std::complex<double> >
{
    %pythonappend CLASS
    %{ 
        self._basisFunctionType = "complex128"
        self._resultType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT

%define BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS);
BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS);
%enddef // BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT

%pythoncode %{

def constructObjectTemplatedOnBasisAndResult(className, basisFunctionType, resultType,
                                             *args, **kwargs):
    fullName = className + "_" + basisFunctionType + "_" + resultType
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
