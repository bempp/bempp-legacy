// Copyright (C) 2011-2015 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef modified_helmholtz_operators_imp_hpp
#define modified_helmholtz_operators_imp_hpp

#include "modified_helmholtz_operators.hpp"

#include "../assembly/blas_quadrature_helper.hpp"
#include "../assembly/context.hpp"
#include "../assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "../assembly/potential_operator.hpp"
#include "../assembly/assembled_potential_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_single_layer_potential_kernel_interpolated_functor.hpp"
#include "../fiber/modified_helmholtz_3d_double_layer_potential_kernel_interpolated_functor.hpp"
#include "../fiber/modified_helmholtz_3d_double_layer_potential_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_adjoint_double_layer_potential_kernel_interpolated_functor.hpp"
#include "../fiber/modified_helmholtz_3d_adjoint_double_layer_potential_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_kernel_interpolated_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_transformation_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_transformation_functor_2.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_integrand_functor_2.hpp"
#include "modified_helmholtz_3d_single_layer_potential_operator.hpp"
#include "modified_helmholtz_3d_double_layer_potential_operator.hpp"

#include "../grid/max_distance.hpp"

#include "../fiber/typical_test_scalar_kernel_trial_integral.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
modifiedHelmholtzSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber,
    const std::string &label, int symmetry)
{
    Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<
      KernelType> NoninterpolatedKernelFunctor;
  typedef Fiber::
      ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<
          KernelType> InterpolatedKernelFunctor;
  typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
      TransformationFunctor;
  typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctorExt<
      BasisFunctionType, KernelType, ResultType, 1> IntegrandFunctor;

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType,
                                                    KernelType, ResultType> Op;
  shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                            ResultType>> integral;
  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange))
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
  else
    integral.reset(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
        IntegrandFunctor()));

  shared_ptr<Op> newOp;
  if (useInterpolation)
    newOp.reset(
        new Op(domain, range, dualToRange, label, symmetry,
               InterpolatedKernelFunctor(
                   waveNumber,
                   1.1 * maxDistance(*domain->grid(), *dualToRange->grid()),
                   interpPtsPerWavelength),
               TransformationFunctor(), TransformationFunctor(), integral));
  else
    newOp.reset(new Op(domain, range, dualToRange, label, symmetry,
                       NoninterpolatedKernelFunctor(waveNumber),
                       TransformationFunctor(), TransformationFunctor(),
                       integral));
  return newOp;

}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
modifiedHelmholtzDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber,
    const std::string &label, int symmetry)
{
  Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::ModifiedHelmholtz3dDoubleLayerPotentialKernelFunctor<
      KernelType> NoninterpolatedKernelFunctor;
  typedef Fiber::
      ModifiedHelmholtz3dDoubleLayerPotentialKernelInterpolatedFunctor<
          KernelType> InterpolatedKernelFunctor;
  typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
      TransformationFunctor;
  typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctorExt<
      BasisFunctionType, KernelType, ResultType, 1> IntegrandFunctor;

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType,
                                                    KernelType, ResultType> Op;

  shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                            ResultType>> integral;
  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange))
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
  else
    integral.reset(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
        IntegrandFunctor()));

  shared_ptr<Op> newOp;
  if (useInterpolation)
    newOp.reset(
        new Op(domain, range, dualToRange, label, symmetry,
               InterpolatedKernelFunctor(
                   waveNumber,
                   1.1 * maxDistance(*domain->grid(), *dualToRange->grid()),
                   interpPtsPerWavelength),
               TransformationFunctor(), TransformationFunctor(), integral));
  else
    newOp.reset(new Op(domain, range, dualToRange, label, symmetry,
                       NoninterpolatedKernelFunctor(waveNumber),
                       TransformationFunctor(), TransformationFunctor(),
                       integral));

  return newOp;

}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
modifiedHelmholtzAdjointDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber,
    const std::string &label, int symmetry)
{

  Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::ModifiedHelmholtz3dAdjointDoubleLayerPotentialKernelFunctor<
      KernelType> NoninterpolatedKernelFunctor;
  typedef Fiber::
      ModifiedHelmholtz3dAdjointDoubleLayerPotentialKernelInterpolatedFunctor<
          KernelType> InterpolatedKernelFunctor;
  typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
      TransformationFunctor;
  typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctorExt<
      BasisFunctionType, KernelType, ResultType, 1> IntegrandFunctor;

  shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                            ResultType>> integral;
  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange))
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
  else
    integral.reset(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
        IntegrandFunctor()));

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType,
                                                    KernelType, ResultType> Op;
  shared_ptr<Op> newOp;
  if (useInterpolation)
    newOp.reset(
        new Op(domain, range, dualToRange, label, symmetry,
               InterpolatedKernelFunctor(
                   waveNumber,
                   1.1 * maxDistance(*domain->grid(), *dualToRange->grid()),
                   interpPtsPerWavelength),
               TransformationFunctor(), TransformationFunctor(), integral));
  else
    newOp.reset(new Op(domain, range, dualToRange, label, symmetry,
                       NoninterpolatedKernelFunctor(waveNumber),
                       TransformationFunctor(), TransformationFunctor(),
                       integral));
  return newOp;

}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
modifiedHelmholtzHypersingularBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber,
    const std::string &label, int symmetry)
{

  Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  int interpPtsPerWavelength = parameterList.get<int>("options.assembly.interpolationPointsPerWavelength");
  bool useInterpolation = parameterList.get<bool>("options.assembly.enableInterpolationForOscillatoryKernels");

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::ModifiedHelmholtz3dHypersingularKernelFunctor<KernelType>
      NoninterpolatedKernelFunctor;
  typedef Fiber::ModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor<
      KernelType> InterpolatedKernelFunctor;
  typedef Fiber::ModifiedHelmholtz3dHypersingularTransformationFunctor<
      CoordinateType> TransformationFunctor;
  typedef Fiber::ModifiedHelmholtz3dHypersingularIntegrandFunctor2<
      BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

  typedef Fiber::ModifiedHelmholtz3dHypersingularTransformationFunctor2<
      CoordinateType> TransformationFunctorWithBlas;

  CoordinateType maxDistance_ =
      static_cast<CoordinateType>(1.1) *
      maxDistance(*domain->grid(), *dualToRange->grid());

  shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                            ResultType>> integral;
  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange)) {
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
  } else {
    integral.reset(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
        IntegrandFunctor()));
  }

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType, KernelType,
                                               ResultType> Op;
  shared_ptr<Op> newOp;
  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange)) {
    shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                              ResultType>> integral;
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
    if (useInterpolation)
      newOp.reset(new Op(
          domain, range, dualToRange, label, symmetry,
          InterpolatedKernelFunctor(waveNumber, maxDistance_,
                                    interpPtsPerWavelength),
          TransformationFunctorWithBlas(), TransformationFunctorWithBlas(),
          integral));     
    else
      newOp.reset(new Op(
          domain, range, dualToRange, label, symmetry,
          NoninterpolatedKernelFunctor(waveNumber),
          TransformationFunctorWithBlas(), TransformationFunctorWithBlas(),
          integral));
  } else { // no blas
    if (useInterpolation)
      newOp.reset(new Op(
          domain, range, dualToRange, label, symmetry,
          InterpolatedKernelFunctor(waveNumber, maxDistance_,
                                    interpPtsPerWavelength),
          TransformationFunctor(), TransformationFunctor(), IntegrandFunctor()));
    else
      newOp.reset(new Op(
          domain, range, dualToRange, label, symmetry,
          NoninterpolatedKernelFunctor(waveNumber), TransformationFunctor(),
          TransformationFunctor(), IntegrandFunctor()));
  }

  return newOp;

}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>> modifiedHelmholtzSingleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        ResultType waveNumber,
        const ParameterList& parameterList)
{

        typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
        shared_ptr<Matrix<CoordinateType>> pointsPtr(
                new Matrix<CoordinateType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, ResultType>> op(
                new ModifiedHelmholtz3dSingleLayerPotentialOperator<BasisFunctionType>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>> modifiedHelmholtzDoubleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        ResultType waveNumber,
        const ParameterList& parameterList)
{

        typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
        shared_ptr<Matrix<CoordinateType>> pointsPtr(
                new Matrix<CoordinateType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, ResultType>> op(
                new ModifiedHelmholtz3dDoubleLayerPotentialOperator<BasisFunctionType>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

}

#endif

