// Macros for classes templated on basis function type and result type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)
    template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
    class CLASS;
%enddef

// Invoke this macro for all *base* classes templated on BasisFunctionType, 
// ResultType and GeometryFactory
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)
    %extend CLASS<float, float, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string resultType() const { return "float32"; }
    }

    %extend CLASS<float, std::complex<float>, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<float>, std::complex<float>, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "complex64"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<double, double, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "float64"; }
    }

    %extend CLASS<double, std::complex<double>, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }

    %extend CLASS<std::complex<double>, std::complex<double>, Bempp::GeometryFactory>
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY

%define BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)
    %template()
        CLASS<float, float, Bempp::GeometryFactory>;
    %template()
        CLASS<float, std::complex<float>, Bempp::GeometryFactory>;
    %template()
        CLASS<std::complex<float>, std::complex<float>, Bempp::GeometryFactory>;

    %template()
        CLASS<double, double, Bempp::GeometryFactory>;
    %template()
        CLASS<double, std::complex<double>, Bempp::GeometryFactory>;
    %template()
        CLASS<std::complex<double>, std::complex<double>, Bempp::GeometryFactory>
%enddef // BEMPP_INSTANTIATE_ANONYMOUSLY_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)
    %template(CLASS ## _float32_float32)
        CLASS<float, float, Bempp::GeometryFactory>;
    %template(CLASS ## _float32_complex64)
        CLASS<float, std::complex<float>, Bempp::GeometryFactory>;
    %template(CLASS ## _complex64_complex64)
        CLASS<std::complex<float>, std::complex<float>, Bempp::GeometryFactory>;

    %template(CLASS ## _float64_float64)
        CLASS<double, double, Bempp::GeometryFactory>;
    %template(CLASS ## _float64_complex128)
        CLASS<double, std::complex<double>, Bempp::GeometryFactory>;
    %template(CLASS ## _complex128_complex128)
        CLASS<std::complex<double>, std::complex<double>, Bempp::GeometryFactory>
%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY
