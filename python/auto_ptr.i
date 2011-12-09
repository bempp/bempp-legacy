#ifdef SWIGPYTHON
%define AUTO_PTR_TYPEMAPS(TYPE...)
%typemap (out) std::auto_ptr<TYPE > %{
    %set_output(SWIG_NewPointerObj(%as_voidptr($1.release()), 
        $descriptor(TYPE *), SWIG_POINTER_OWN | %newpointer_flags));
%}
%template() std::auto_ptr<TYPE >;
%enddef
#endif

namespace std {
   template <class T> class auto_ptr {};
}

/* AUTO_PTR_TYPEMAPS(MyStruct) */

/* %inline %{ */
/* #include <memory> */
/* struct MyStruct {}; */
/* std::auto_ptr<MyStruct> factory_function() { */
/*    return std::auto_ptr<MyStruct>(new MyStruct()); */
/* } */
/* %} */
