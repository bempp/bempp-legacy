// Macros for classes templated on basis function, kernel and result type

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(CLASS)
    %template(CLASS ## _float32_float32_float32)
        CLASS<float, float, float>;
    %template(CLASS ## _float32_float32_complex64)
        CLASS<float, float, std::complex<float> >;
    %template(CLASS ## _float32_complex64_complex64)
        CLASS<float, std::complex<float>, std::complex<float> >;
    %template(CLASS ## _complex64_float32_complex64)
        CLASS<std::complex<float>, float, std::complex<float> >;
    %template(CLASS ## _complex64_complex64_complex64)
        CLASS<std::complex<float>, std::complex<float>, std::complex<float> >;

    %template(CLASS ## _float64_float64_float64)
        CLASS<double, double, double>;
    %template(CLASS ## _float64_float64_complex128)
        CLASS<double, double, std::complex<double> >;
    %template(CLASS ## _float64_complex128_complex128)
        CLASS<double, std::complex<double>, std::complex<double> >;
    %template(CLASS ## _complex128_float64_complex128)
        CLASS<std::complex<double>, double, std::complex<double> >;
    %template(CLASS ## _complex128_complex128_complex128)
        CLASS<std::complex<double>, std::complex<double>, std::complex<double> >
%enddef
