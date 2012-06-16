// Macros for classes templated on basis function type

%define BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(CLASS)
template <typename BasisFunctionType> class CLASS;
%enddef

%define BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(CLASS)

%template(CLASS ## _float32)
    CLASS<float>;
%template(CLASS ## _complex64)
    CLASS<std::complex<float> >;
%template(CLASS ## _float64)
    CLASS<double>;
%template(CLASS ## _complex128)
    CLASS<std::complex<double> >;

%enddef // BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS

// deprecated
%define BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS(CLASS)

%extend CLASS<float>
{
    %pythonappend CLASS
    %{
        self._basisFunctionType = "float32"
    %}
}

%extend CLASS<double>
{
    %pythonappend CLASS
    %{
        self._basisFunctionType = "float64"
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythonappend CLASS
    %{
        self._basisFunctionType = "complex64"
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythonappend CLASS
    %{
        self._basisFunctionType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS

%define BEMPP_PYTHON_EXTEND_INTERFACE_CLASS_TEMPLATED_ON_BASIS(CLASS)

%extend CLASS<float>
{
    %pythoncode
    %{
        def _initTypes(self):
            self._basisFunctionType = "float32"
    %}
}

%extend CLASS<double>
{
    %pythoncode
    %{
        def _initTypes(self):
            self._basisFunctionType = "float64"
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythoncode
    %{
        def _initTypes(self):
            self._basisFunctionType = "complex64"
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythoncode
    %{
        def _initTypes(self):
            self._basisFunctionType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_INTERFACE_CLASS_TEMPLATED_ON_BASIS

%define BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(CLASS)

%extend CLASS<float>
{
    %pythonappend CLASS
    %{
        self._initTypes()
    %}
}

%extend CLASS<double>
{
    %pythonappend CLASS
    %{
        self._initTypes()
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythonappend CLASS
    %{
        self._initTypes()
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythonappend CLASS
    %{
        self._initTypes()
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS

%define BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS(CLASS)
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(CLASS);
BEMPP_PYTHON_EXTEND_CLASS_TEMPLATED_ON_BASIS(CLASS);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(CLASS);
%enddef // BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS

%define BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS(CLASS)
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(CLASS);
BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(CLASS);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(CLASS);
%enddef // BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS

%pythoncode %{

def constructObjectTemplatedOnBasis(className, basisFunctionType,*args, **kwargs):
    fullName = className + "_" + checkType(basisFunctionType)
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
