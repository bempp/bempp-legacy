// Automatically specify argument lists and return types in docstrings
%feature("autodoc", 0);

%define DECLARE_CLASS_DOCSTRING(class)
%feature("docstring", class ## _docstring) class
%enddef

// class -> class name
// method -> method name 
// auto_param_list -> if true, the parameter list and result's type will be 
// generated automatically, otherwise it should be provided explicitly 
// in the docstring
%define DECLARE_METHOD_DOCSTRING(class, method, auto_param_list)
#if auto_param_list
    %feature("autodoc", 0) class::method;
#endif
%feature("docstring", class ## _ ## method ## _docstring) class::method;
%enddef

%define DECLARE_TEMPLATE_CLASS_DOCSTRING(python_class, cpp_class)
%feature("docstring", python_class ## _docstring) cpp_class;
%enddef

%define DECLARE_TEMPLATE_METHOD_DOCSTRING(python_class, cpp_class, method, auto_param_list)
#if !auto_param_list
    %feature("autodoc", python_class ## _ ## method ## _autodoc_docstring) cpp_class::method;
#endif
%feature("docstring", python_class ## _ ## method ## _docstring) cpp_class::method;
%enddef

%define DECLARE_TEMPLATE_VALUE_METHOD_AUTO_DOCSTRING(class, method)
%feature("autodoc", 2) class<float>::method;
%feature("autodoc", 2) class<double>::method;
%feature("autodoc", 2) class<std::complex<float> >::method;
%feature("autodoc", 2) class<std::complex<double> >::method;
%enddef

%define BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS(class)
%feature("docstring", class ## _docstring(float32)) 
         Bempp:: ## class ## < float >;
%feature("docstring", class ## _docstring(complex64))
         Bempp:: ## class ## < std::complex<float> >;
%feature("docstring", class ## _docstring(float64))
         Bempp:: ## class ## < double >;
%feature("docstring", class ## _docstring(complex128))
         Bempp:: ## class ## < std::complex<double> >;
%enddef

%define BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS(class, method)
%feature("docstring", class ## _ ## method ## _docstring(float32)) 
         Bempp:: class < float > :: method;
%feature("docstring", class ## _ ## method ## _docstring(complex64))
         Bempp:: class < std::complex<float> > :: method;
%feature("docstring", class ## _ ## method ## _docstring(float64))
         Bempp:: class < double > :: method;
%feature("docstring", class ## _ ## method ## _docstring(complex128))
         Bempp:: class < std::complex<double> > :: method;
%enddef

%define BEMPP_DECLARE_DOCSTRING_FOR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(class)
%feature("docstring", class ## _docstring(float32, float32)) 
         Bempp:: ## class ## < float, float >;
%feature("docstring", class ## _docstring(float32, complex64))
         Bempp:: ## class ## < float, std::complex<float> >;
%feature("docstring", class ## _docstring(complex64, complex64))
         Bempp:: ## class ## < std::complex<float>, std::complex<float> >;
%feature("docstring", class ## _docstring(float64, float64))
         Bempp:: ## class ## < double, double >;
%feature("docstring", class ## _docstring(float64, complex128))
         Bempp:: ## class ## < double, std::complex<double> >;
%feature("docstring", class ## _docstring(complex128, complex128))
         Bempp:: ## class ## < std::complex<double>, std::complex<double> >;
%enddef

%define BEMPP_DECLARE_DOCSTRING_FOR_METHOD_OF_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(class, method)
%feature("docstring", class ## _ ## method ## _docstring(float32, float32)) 
         Bempp:: class < float, float > :: method;
%feature("docstring", class ## _ ## method ## _docstring(float32, complex64))
         Bempp:: class < float, std::complex<float> > :: method;
%feature("docstring", class ## _ ## method ## _docstring(complex64, complex64))
         Bempp:: class < std::complex<float>, std::complex<float> > :: method;
%feature("docstring", class ## _ ## method ## _docstring(float64, float64))
         Bempp:: class < double, double > :: method;
%feature("docstring", class ## _ ## method ## _docstring(float64, complex128))
         Bempp:: class < double, std::complex<double> > :: method;
%feature("docstring", class ## _ ## method ## _docstring(complex128, complex128))
         Bempp:: class < std::complex<double>, std::complex<double> > :: method;
%enddef
