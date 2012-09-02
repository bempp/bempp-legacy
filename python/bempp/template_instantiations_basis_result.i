// Macros for classes templated on basis function type and result type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
    template <typename BasisFunctionType, typename ResultType> class CLASS;
%enddef

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
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
%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT

// Invoke this macro for all *base* classes templated on BasisFunctionType and
// ResultType
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
    %extend CLASS<float, float>
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string resultType() const { return "float32"; }
    }

    %extend CLASS<float, std::complex<float> >
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<float>, std::complex<float> >
    {
        std::string basisFunctionType() const { return "complex64"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<double, double>
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "float64"; }
    }

    %extend CLASS<double, std::complex<double> >
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }

    %extend CLASS<std::complex<double>, std::complex<double> >
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT

%define BEMPP_DECLARE_SHARED_PTR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CLASS)
    %shared_ptr(CLASS<float, float>);
    %shared_ptr(CLASS<float, std::complex<float> >);
    %shared_ptr(CLASS<std::complex<float>, std::complex<float> >);
    %shared_ptr(CLASS<double, double>);
    %shared_ptr(CLASS<double, std::complex<double> >)
    %shared_ptr(CLASS<std::complex<double>, std::complex<double> >)
%enddef            
