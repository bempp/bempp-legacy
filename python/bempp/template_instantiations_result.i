// Macros for classes templated on result type

%define BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_RESULT(CLASS)
template <typename ResultType> class CLASS;
%enddef

%define BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_RESULT(CLASS)
    %template(CLASS ## _float32)
        CLASS<float>;
    %template(CLASS ## _complex64)
        CLASS<std::complex<float> >;

    %template(CLASS ## _float64)
        CLASS<double>;
    %template(CLASS ## _complex128)
        CLASS<std::complex<double> >;
%enddef // BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_RESULT

// Invoke this macro for all *base* classes templated on ResultType
%define BEMPP_EXTEND_CLASS_TEMPLATED_ON_RESULT(CLASS)
    %extend CLASS<float>
    {
        std::string resultType() const { return "float32"; }
    }

    %extend CLASS<double>
    {
        std::string resultType() const { return "float64"; }
    }

    %extend CLASS<std::complex<float> >
    {
        std::string resultType() const { return "complex64"; }
    }

    %extend CLASS<std::complex<double> >
    {
        std::string resultType() const { return "complex128"; }
    }
%enddef // BEMPP_EXTEND_CLASS_TEMPLATED_ON_RESULT

%pythoncode 
%{

def constructObjectTemplatedOnResult(className, resultType,
                                     *args, **kwargs):
    fullName = className + "_" + resultType
    try:
        class_ = globals()[fullName]
    except KeyError:
        raise TypeError("Class " + fullName + " does not exist.")
    return class_(*args, **kwargs)

%}
