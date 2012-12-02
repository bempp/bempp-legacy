%{
#include "assembly/l2_norm.hpp"
%}

%inline %{

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType
L2NormOfDifferenceFromPythonSurfaceNormalIndependentFunctor(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const PythonSurfaceNormalIndependentFunctor<ResultType>& functor,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    return L2NormOfDifference(gridFunction,
                              surfaceNormalIndependentFunction(functor),
                              quadStrategy,
                              options);
}

template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType
L2NormOfDifferenceFromPythonSurfaceNormalDependentFunctor(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const PythonSurfaceNormalDependentFunctor<ResultType>& functor,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
{
    return L2NormOfDifference(gridFunction,
                              surfaceNormalDependentFunction(functor),
                              quadStrategy,
                              options);
}

} // namespace Bempp

%}

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(L2NormOfDifferenceFromPythonSurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(L2NormOfDifferenceFromPythonSurfaceNormalDependentFunctor);

} // namespace Bempp