%define UNIQUE_PTR_TYPEMAPS(TYPE...)
%typemap (out) std::unique_ptr<TYPE > %{
    %set_output(SWIG_NewPointerObj(%as_voidptr($1.release()), 
        $descriptor(TYPE *), SWIG_POINTER_OWN | %newpointer_flags));
%}
%template() std::unique_ptr<TYPE >;
%enddef

%define UNIQUE_PTR_TYPEMAPS_FOR_CLASS_TEMPLATED_ON_RESULT(TYPE)
UNIQUE_PTR_TYPEMAPS(TYPE<float>);
UNIQUE_PTR_TYPEMAPS(TYPE<double>);
UNIQUE_PTR_TYPEMAPS(TYPE<std::complex<float> >);
UNIQUE_PTR_TYPEMAPS(TYPE<std::complex<double> >);
%enddef

%define UNIQUE_PTR_TYPEMAPS_FOR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(TYPE)
UNIQUE_PTR_TYPEMAPS(TYPE<float, float>);
UNIQUE_PTR_TYPEMAPS(TYPE<float, std::complex<float> >);
UNIQUE_PTR_TYPEMAPS(TYPE<std::complex<float>, std::complex<float> >);
UNIQUE_PTR_TYPEMAPS(TYPE<double, double>);
UNIQUE_PTR_TYPEMAPS(TYPE<double, std::complex<double> >);
UNIQUE_PTR_TYPEMAPS(TYPE<std::complex<double>, std::complex<double> >);
%enddef

namespace std {
   template <class T> class unique_ptr {};
}
