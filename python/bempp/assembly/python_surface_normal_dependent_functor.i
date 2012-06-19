%{
#include "assembly/surface_normal_dependent_functor.hpp"
#include <armadillo>
  %}


%typemap(in) PyObject *pyfunc {
  if (!PyCallable_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "Callable object expected!");
    return NULL;
  }
  $1=$input;
 }


%feature("pythonappend") Bempp::PythonSurfaceNormalDependentFunctor_float32::PythonSurfaceNormalDependentFunctor_float32(PyObject *pyFunc,int argumentDimension, int resultDimension)
 %{
self._valueType="float32"
  %}
%feature("pythonappend") Bempp::PythonSurfaceNormalDependentFunctor_float64::PythonSurfaceNormalDependentFunctor_float64(PyObject *pyFunc,int argumentDimension, int resultDimension)
 %{
self._valueType="float64"
  %}
%feature("pythonappend") Bempp::PythonSurfaceNormalDependentFunctor_complex64::PythonSurfaceNormalDependentFunctor_complex64(PyObject *pyFunc,int argumentDimension, int resultDimension)
 %{
self._valueType="complex64"
  %}
%feature("pythonappend") Bempp::PythonSurfaceNormalDependentFunctor_complex128::PythonSurfaceNormalDependentFunctor_complex128(PyObject *pyFunc,int argumentDimension, int resultDimension)
 %{
self._valueType="complex128"
  %}



%define BEMPP_DECLARE_PYTHON_SURFACE_NORMAL_DEPENDENT_FUNCTOR( TYPE , NPY_NAME , NPY_TYPE )



%inline %{

namespace Bempp{
class PythonSurfaceNormalDependentFunctor_## NPY_NAME : public SurfaceNormalDependentFunctor< TYPE >
{
 public:
  PythonSurfaceNormalDependentFunctor_## NPY_NAME(PyObject *pyFunc,int argumentDimension, int resultDimension) :
  m_pyFunc(pyFunc), m_argumentDimension(argumentDimension), m_resultDimension(resultDimension)
{
  if (!PyCallable_Check(pyFunc)) PyErr_SetString(PyExc_TypeError, "Python object is not callable");    
  Py_INCREF(m_pyFunc); // Increase shared pointer
}

  int argumentDimension() const {
    return m_argumentDimension;
  }

  int resultDimension() const {
    return m_resultDimension;
  }

  void evaluate(const arma::Col<CoordinateType>& point,
                const arma::Col<CoordinateType>& normal,		
		arma::Col<ValueType>& result_) const
  {
    // Create the input array    
    npy_intp dims1[1];
    dims1[0]=point.n_rows;
    PyObject* pyPoint = PyArray_ZEROS(1, dims1, NPY_TYPE , NPY_FORTRAN);
    TYPE* pdata=(TYPE*)array_data(pyPoint);
    for (size_t i=0;i<dims1[0];i++) pdata[i]=point(i);

    npy_intp dims2[1];
    dims2[0]=normal.n_rows;
    PyObject* pyNormal = PyArray_ZEROS(1, dims2, NPY_TYPE , NPY_FORTRAN);
    TYPE* ndata=(TYPE*)array_data(pyNormal);
    for (size_t i=0;i<dims2[0];i++) ndata[i]=normal(i);
    
    // Call into Python
    PyObject* pyReturnVal= PyObject_CallFunctionObjArgs(m_pyFunc,pyPoint, pyNormal, NULL);
    if (!pyReturnVal) {
      Py_XDECREF(pyPoint);
      Py_XDECREF(pyNormal);
      throw std::runtime_error("Callable did not execute successfully");
    }
    PyObject* pyReturnValArray=PyArray_FROM_OT(pyReturnVal,NPY_TYPE);
    if (!pyReturnValArray) {
      Py_XDECREF(pyPoint);
      Py_XDECREF(pyNormal);
      Py_XDECREF(pyReturnVal);
      throw std::runtime_error("Result from callable cannot be converted to array.");
    }
    int is_new_object;
    PyArrayObject* pyReturnValArrayCont=obj_to_array_contiguous_allow_conversion(pyReturnValArray,NPY_TYPE,&is_new_object);

    // Check size of array
    int asize;
    if (PyArray_CheckScalar(pyReturnValArrayCont)){
      asize=1;
    }
    else{
      // Check number of dimensions
      if (array_numdims(pyReturnValArrayCont)!=1){
	Py_XDECREF(pyPoint);
	Py_XDECREF(pyNormal);
	Py_XDECREF(pyReturnVal);
	Py_XDECREF(pyReturnValArray);
	Py_XDECREF(pyReturnValArrayCont);
	throw std::runtime_error("Return array has wrong dimensions!");
      }
      asize=array_size(pyReturnValArrayCont,0);
    }
    if (asize!=m_resultDimension) {
      Py_XDECREF(pyNormal);
      Py_XDECREF(pyPoint);
      Py_XDECREF(pyReturnVal);
      Py_XDECREF(pyReturnValArray);
      Py_XDECREF(pyReturnValArrayCont);
      throw std::runtime_error("Return array has wrong dimensions");
    }
     
    // Copy data back
    TYPE* data= (TYPE*) array_data(pyReturnValArrayCont);
    for (size_t i=0;i<m_resultDimension;i++) result_(i)=data[i];

    // Clean Up

    Py_XDECREF(pyPoint);
    Py_XDECREF(pyNormal);
    Py_XDECREF(pyReturnVal);
    Py_XDECREF(pyReturnValArray);
    Py_XDECREF(pyReturnValArrayCont);

  }    


   ~PythonSurfaceNormalDependentFunctor_##NPY_NAME(){
    Py_DECREF(m_pyFunc);
  }
 private:
    PyObject* m_pyFunc;
    int m_argumentDimension;
    int m_resultDimension;
};
}


%}

%enddef

BEMPP_DECLARE_PYTHON_SURFACE_NORMAL_DEPENDENT_FUNCTOR(float,float32,NPY_FLOAT)
BEMPP_DECLARE_PYTHON_SURFACE_NORMAL_DEPENDENT_FUNCTOR(double,float64,NPY_DOUBLE)
BEMPP_DECLARE_PYTHON_SURFACE_NORMAL_DEPENDENT_FUNCTOR(std::complex<float>,complex64,NPY_CFLOAT)
BEMPP_DECLARE_PYTHON_SURFACE_NORMAL_DEPENDENT_FUNCTOR(std::complex<double>,complex128,NPY_CDOUBLE)

  %pythoncode %{

def surfaceNormalDependentFunctor(fun,valueType='float64',argumentDimension=3,resultDimension=1):

    dtype=checkType(valueType)
    base=None
    if dtype=='float32':
        functor=PythonSurfaceNormalDependentFunctor_float32
    elif dtype=='float64':
        functor=PythonSurfaceNormalDependentFunctor_float64
    elif dtype=='complex64':
        functor=PythonSurfaceNormalDependentFunctor_complex64
    elif dtype=='complex128':
        functor=PythonSurfaceNormalDependentFunctor_complex128
    return functor(fun,argumentDimension,resultDimension);

	      %}

