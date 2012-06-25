// Macros for classes templated on basis function type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(CLASS)
    template <typename BasisFunctionType> class CLASS;
%enddef

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(CLASS)
    %template(CLASS ## _float32)
        CLASS<float>;
    %template(CLASS ## _complex64)
        CLASS<std::complex<float> >;
    %template(CLASS ## _float64)
        CLASS<double>;
    %template(CLASS ## _complex128)
        CLASS<std::complex<double> >;
%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS

// Invoke this macro for all *base* classes templated on BasisFunctionType
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS(CLASS)
    %extend CLASS<float>
    {
        std::string basisFunctionType() const { return "float32"; }
    }

    %extend CLASS<double>
    {
        std::string basisFunctionType() const { return "float64"; }
    }

    %extend CLASS<std::complex<float> >
    {
        std::string basisFunctionType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<double> >
    {
        std::string basisFunctionType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS

%pythoncode %{

def constructObjectTemplatedOnBasis(className, basisFunctionType, *args, **kwargs):
    fullName = className + "_" + checkType(basisFunctionType)
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
