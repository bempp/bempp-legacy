#ifndef PY_FUNCTORS_HPP
#define PY_FUNCTORS_HPP

#include "bempp/fiber/surface_normal_dependent_function.hpp"
#include "bempp/fiber/scalar_traits.hpp"
#include <armadillo>
namespace Bempp
{

template <typename ValueType_>
class PythonSurfaceNormalDependentFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef void (*pyFunc_t)(const arma::Col<CoordinateType>&, const arma::Col<CoordinateType>&,
            arma::Col<ValueType>&);

    PythonSurfaceNormalDependentFunctor(
        pyFunc_t pyFunc,
        int argumentDimension, int resultDimension) :
            m_pyFunc(pyFunc),
            m_argumentDimension(argumentDimension),
            m_resultDimension(resultDimension){} 

    int argumentDimension() const {
        return m_argumentDimension;
    }

    int resultDimension() const {
        return m_resultDimension;
    }

    void evaluate(const arma::Col<CoordinateType>& point, const arma::Col<CoordinateType>& normal,
                  arma::Col<ValueType>& result_) const
    {
        m_pyFunc(point,normal,result_);
    }

private:
    pyFunc_t m_pyFunc;
    int m_argumentDimension;
    int m_resultDimension;
};


template <typename ValueType>
shared_ptr<Fiber::Function<ValueType>> _py_surface_normal_dependent_function(
        typename PythonSurfaceNormalDependentFunctor<ValueType>::pyFunc_t pyFunc, 
        int argumentDimension, int resultDimension)
{
    return shared_ptr<Fiber::Function<ValueType>>(
        new Fiber::SurfaceNormalDependentFunction<PythonSurfaceNormalDependentFunctor<ValueType>>(
            PythonSurfaceNormalDependentFunctor<ValueType>(pyFunc,argumentDimension,resultDimension)));
}
} // namespace Bempp


#endif
