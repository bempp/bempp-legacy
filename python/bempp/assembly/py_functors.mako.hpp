#ifndef PY_FUNCTORS_HPP
#define PY_FUNCTORS_HPP

#include "bempp/fiber/surface_normal_and_domain_index_dependent_function.hpp"
#include "bempp/fiber/scalar_traits.hpp"
#include <vector>
#include <armadillo>
#include <Python.h>
#include <numpy/arrayobject.h>

namespace Bempp
{

    template <typename T> struct NumpyType {
    };

    template <> struct NumpyType<float> {
        enum { value = 11 };
    };

    template <> struct NumpyType<double> {
        enum { value = 12 };
    };

    template <> struct NumpyType<std::complex<float>> {
        enum { value = 14 };
    };

    template <> struct NumpyType<std::complex<double>> {
        enum { value = 15 };
    };


template <typename ValueType_>
class PythonFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef void (*pyFunc_t)(PyObject* x, PyObject* normal, int domainIndex, PyObject* result, PyObject* callable);



    PythonFunctor(
        pyFunc_t pyFunc, PyObject* callable,
        int argumentDimension, int resultDimension) :
            m_pyFunc(pyFunc),
            m_callable(callable),
            m_argumentDimension(argumentDimension),
            m_resultDimension(resultDimension)
            {

            Py_INCREF(m_callable);

            npy_intp pyArgumentDimension = m_argumentDimension;
            npy_intp pyResultDimension = m_resultDimension;

            m_x = PyArray_ZEROS(1,&pyArgumentDimension,NumpyType<CoordinateType>::value,1);
            m_normal = PyArray_ZEROS(1,&pyArgumentDimension,NumpyType<CoordinateType>::value,1);
            m_result = PyArray_ZEROS(1,&pyResultDimension,NumpyType<ValueType>::value,1);

            } 

    PythonFunctor(const PythonFunctor<ValueType>& other):
        m_pyFunc(other.m_pyFunc), m_argumentDimension(other.m_argumentDimension),
        m_resultDimension(other.m_resultDimension),m_x(other.m_x),m_normal(other.m_normal),
        m_result(other.m_result),m_callable(other.m_callable) {

            Py_INCREF(m_callable);
            Py_INCREF(m_x);
            Py_INCREF(m_normal);
            Py_INCREF(m_result);
            

        }

    ~PythonFunctor(){

        Py_DECREF(m_x);
        Py_DECREF(m_normal);
        Py_DECREF(m_result);
        Py_DECREF(m_callable);

    }


    int argumentDimension() const {
        return m_argumentDimension;
    }

    int resultDimension() const {
        return m_resultDimension;
    }

    void evaluate(const arma::Col<CoordinateType>& point, const arma::Col<CoordinateType>& normal,
                  int domainIndex, arma::Col<ValueType>& result_) const
    {

        CoordinateType* xPtr = (CoordinateType*)PyArray_DATA(m_x);
        for (int i = 0; i< m_argumentDimension;++i) xPtr[i] = point.at(i);

        CoordinateType* normalPtr = (CoordinateType*)PyArray_DATA(m_normal);
        for (int i = 0; i< m_argumentDimension;++i) normalPtr[i] = normal.at(i);

        m_pyFunc(m_x,m_normal,domainIndex,m_result,m_callable);

        ValueType* resPtr = (ValueType*)PyArray_DATA(m_result);
        for (int i = 0; i< m_resultDimension;++i) result_.at(i) = resPtr[i];
    }

private:
    pyFunc_t m_pyFunc;
    int m_argumentDimension;
    int m_resultDimension;
    PyObject* m_x;
    PyObject* m_normal;
    PyObject* m_result;
    PyObject* m_callable;

};


template <typename ValueType>
shared_ptr<Fiber::Function<ValueType>> _py_surface_normal_dependent_function(
        typename PythonFunctor<ValueType>::pyFunc_t pyFunc,PyObject* callable, 
        int argumentDimension, int resultDimension)
{
    return shared_ptr<Fiber::Function<ValueType>>(
        new Fiber::SurfaceNormalAndDomainIndexDependentFunction<PythonFunctor<ValueType>>(
            PythonFunctor<ValueType>(pyFunc,callable,argumentDimension,resultDimension)));
}
} // namespace Bempp


#endif
