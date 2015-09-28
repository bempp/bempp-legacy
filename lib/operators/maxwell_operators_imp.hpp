

#ifndef bempp_maxwell_operators_imp_hpp
#define bempp_maxwell_operators_imp_hpp

#include "maxwell_operators.hpp"

#include "../assembly/blas_quadrature_helper.hpp"
#include "../assembly/context.hpp"
#include "../assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "../assembly/potential_operator.hpp"
#include "../assembly/assembled_potential_operator.hpp"
#include "../assembly/maxwell_3d_double_layer_potential_operator.hpp"
#include "../assembly/maxwell_3d_single_layer_potential_operator.hpp"
#include "../assembly/maxwell_3d_far_field_double_layer_potential_operator.hpp"
#include "../assembly/maxwell_3d_far_field_single_layer_potential_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_kernel_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_kernel_interpolated_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_operators_transformation_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_integrand_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_operators_kernel_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_operators_kernel_interpolated_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_boundary_operator_integrand_functor.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"


#include "../grid/max_distance.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {


template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
maxwellElectricFieldBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string& label, int symmetry)
{

  typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
  typedef typename ScalarTraits<BasisFunctionType>::ComplexType ResultType;
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef Fiber::ModifiedMaxwell3dSingleLayerBoundaryOperatorKernelFunctor<
      KernelType> KernelFunctor;
  typedef Fiber::
      ModifiedMaxwell3dSingleLayerBoundaryOperatorKernelInterpolatedFunctor<
          KernelType> KernelInterpolatedFunctor;
  typedef Fiber::ModifiedMaxwell3dSingleLayerOperatorsTransformationFunctor<
      CoordinateType> TransformationFunctor;
  typedef Fiber::ModifiedMaxwell3dSingleLayerBoundaryOperatorIntegrandFunctor<
      BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType,
                                                    KernelType, ResultType> Op;
  if (useInterpolation)
    return 
        boost::make_shared<Op>(
            domain, range, dualToRange, label, symmetry,
            KernelInterpolatedFunctor(
                waveNumber / KernelType(0., 1.),
                1.1 * maxDistance(*domain->grid(), *dualToRange->grid()),
                interpPtsPerWavelength),
            TransformationFunctor(), TransformationFunctor(),
            IntegrandFunctor());
  else
        return boost::make_shared<Op>(domain, range, dualToRange, label, symmetry,
                               KernelFunctor(waveNumber / KernelType(0., 1.)),
                               TransformationFunctor(), TransformationFunctor(),
                               IntegrandFunctor());
}


template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
maxwellMagneticFieldBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label, int symmetry)
{


  typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
  typedef typename ScalarTraits<BasisFunctionType>::ComplexType ResultType;
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef Fiber::ModifiedMaxwell3dDoubleLayerOperatorsKernelFunctor<KernelType>
      KernelFunctor;
  typedef Fiber::ModifiedMaxwell3dDoubleLayerOperatorsKernelInterpolatedFunctor<
      KernelType> KernelInterpolatedFunctor;
  typedef Fiber::HdivFunctionValueFunctor<CoordinateType> TransformationFunctor;
  typedef Fiber::ModifiedMaxwell3dDoubleLayerBoundaryOperatorIntegrandFunctor<
      BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType,
                                                    KernelType, ResultType> Op;
  if (useInterpolation)
    return 
        boost::make_shared<Op>(
            domain, range, dualToRange, label, symmetry,
            KernelInterpolatedFunctor(
                waveNumber / KernelType(0., 1.),
                1.1 * maxDistance(*domain->grid(), *dualToRange->grid()),
                interpPtsPerWavelength),
            TransformationFunctor(), TransformationFunctor(),
            IntegrandFunctor());
  else
    return 
        boost::make_shared<Op>(domain, range, dualToRange, label, symmetry,
                               KernelFunctor(waveNumber / KernelType(0., 1.)),
                               TransformationFunctor(), TransformationFunctor(),
                               IntegrandFunctor());

}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> 
electricFieldPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{
        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> 
magneticFieldPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<double, std::complex<double>>> op(
                new Maxwell3dDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();


}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> 
electricFieldFarFieldOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> op(
                new Maxwell3dFarFieldSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> 
magneticFieldFarFieldOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<double>> pointsPtr(
                new Matrix<double>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> op(
                new Maxwell3dFarFieldDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

}
#endif
