// Macros for classes templated on value type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(CLASS)
template <typename ValueType> class CLASS;
%enddef

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(CLASS)

%template(CLASS ## _float32)
    CLASS<float>;
%template(CLASS ## _complex64)
    CLASS<std::complex<float> >;

%template(CLASS ## _float64)
    CLASS<double>;
%template(CLASS ## _complex128)
    CLASS<std::complex<double> >;

%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE

// Invoke this macro for all *base* classes templated on ValueType
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(CLASS)
    %extend CLASS<float>
    {
        std::string valueType() const { return "float32"; }
    }

    %extend CLASS<double>
    {
        std::string valueType() const { return "float64"; }
    }

    %extend CLASS<std::complex<float> >
    {
        std::string valueType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<double> >
    {
        std::string valueType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE

%pythoncode 
%{


%}
