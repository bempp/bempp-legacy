// Macros for classes templated on basis function type and result type

%define BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)
template <typename BasisFunctionType, typename ResultType, typename GeometryFactory> class CLASS;
%enddef

%define BEMPP_PYTHON_INSTANTIATE_ANONYMOUSLY_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY(CLASS)

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

%enddef // BEMPP_PYTHON_INSTANTIATE_ANONYMOUSLY_CLASS_TEMPLATED_ON_BASIS_RESULT_AND_GEOMETRY_FACTORY
