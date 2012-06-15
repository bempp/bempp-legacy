%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

//namespace Bempp
//{
//  template <typename BasisFunctionType> class PiecewiseConstantScalarSpace : 
//}


%include "space/piecewise_constant_scalar_space.hpp"

namespace Bempp
{

    BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);

/* template <typename BasisFunctionType> class PiecewiseConstantScalarSpace; */

/* %extend PiecewiseConstantScalarSpace<float>  */
/* { */
/*     %pythonappend PiecewiseConstantScalarSpace  */
/*     %{ */
/*         self._basisFunctionType = "float32" */
/*     %} */
/* } */

/*   %template(PiecewiseConstantScalarSpace_float64)     PiecewiseConstantScalarSpace<double>; */
/*   %template(PiecewiseConstantScalarSpace_float32)     PiecewiseConstantScalarSpace<float>; */
/*   %template(PiecewiseConstantScalarSpace_complex64)   PiecewiseConstantScalarSpace<std::complex<float> >; */
/*   %template(PiecewiseConstantScalarSpace_complex128)  PiecewiseConstantScalarSpace<std::complex<double> >; */

} // namespace Bempp

/* def constructObjectTemplatedOnBasis(className, basisFunctionType, */
/*                                     *args, **kwargs): */
/*     fullName = className + "_" + basisFunctionType */
/*     try: */
/*         print globals().keys() */
/*         class_ = globals()[fullName] */
/*     except KeyError: */
/*         raise TypeError("Class " + fullName + " does not exist.") */
/*     return class_(*args, **kwargs) */

%pythoncode 
%{

def piecewiseConstantScalarSpace(basisFunctionType, grid):
    """Construct a PiecewiseConstantScalarSpace object with a given basisFunctionType"""
    return constructObjectTemplatedOnBasis(
        "PiecewiseConstantScalarSpace", basisFunctionType, 
        grid)

%}

/* %pythoncode %{ */

/* class PiecewiseConstantScalarSpace(Template1,ScalarSpace): */
/*     """Space of piecewise constant scalar functions""" */
/*     def __init__(self,dtype1,*args,**kwargs): */
/*         super(PiecewiseConstantScalarSpace,self).__init__('PiecewiseConstantScalarSpace',dtype1,*args,**kwargs) */

/* 	  %} */

