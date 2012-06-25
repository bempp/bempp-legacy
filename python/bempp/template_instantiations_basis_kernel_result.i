// Macros for classes templated on basis function, kernel and result type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(CLASS)
    template <typename BasisFunctionType, typename KernelType, typename ResultType> class CLASS;
%enddef

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

// Invoke this macro for all *base* classes templated on BasisFunctionType and
// ResultType
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(CLASS)
    %extend CLASS<float, float, float>
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string kernelType() const { return "float32"; }
        std::string resultType() const { return "float32"; }
    }

    %extend CLASS<float, float, std::complex<float> >
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string kernelType() const { return "float32"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<float, std::complex<float>, std::complex<float> >
    {
        std::string basisFunctionType() const { return "float32"; }
        std::string kernelType() const { return "complex64"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<float>, float, std::complex<float> >
    {
        std::string basisFunctionType() const { return "complex64"; }
        std::string kernelType() const { return "float32"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<float>, std::complex<float>, std::complex<float> >
    {
        std::string basisFunctionType() const { return "complex64"; }
        std::string kernelType() const { return "complex64"; }
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<double, double, double>
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string kernelType() const { return "float64"; }
        std::string resultType() const { return "float64"; }
    }

    %extend CLASS<double, double, std::complex<double> >
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string kernelType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }

    %extend CLASS<double, std::complex<double>, std::complex<double> >
    {
        std::string basisFunctionType() const { return "float64"; }
        std::string kernelType() const { return "complex128"; }
        std::string resultType() const { return "complex128"; }
    }

    %extend CLASS<std::complex<double>, double, std::complex<double> >
    {
        std::string basisFunctionType() const { return "complex128"; }
        std::string kernelType() const { return "float64"; }
        std::string resultType() const { return "complex128"; }
    }

    %extend CLASS<std::complex<double>, std::complex<double>, std::complex<double> >
    {
        std::string basisFunctionType() const { return "complex128"; }
        std::string kernelType() const { return "complex128"; }
        std::string resultType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT

%pythoncode %{

# note: None values for types are not supported
def constructObjectTemplatedOnBasisKernelAndResult(
        className, basisFunctionType, kernelType, resultType, *args, **kwargs):

    fullName = (className + "_" +
                checkType(basisFunctionType) + "_" +
                checkType(kernelType) + "_" +
                checkType(resultType))
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
