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

#ifndef laplace_operators_imp_hpp
#define laplace_operators_imp_hpp

#include "laplace_operators.hpp"

#include "../assembly/blas_quadrature_helper.hpp"
#include "../assembly/context.hpp"
#include "../assembly/general_elementary_singular_integral_operator_imp.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"
#include "../assembly/potential_operator.hpp"
#include "../assembly/assembled_potential_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/laplace_3d_adjoint_double_layer_potential_kernel_functor.hpp"
#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/laplace_3d_double_layer_potential_kernel_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"
#include "../fiber/scalar_traits.hpp"

#include "../fiber/typical_test_scalar_kernel_trial_integral.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
{
    Context<BasisFunctionType, ResultType> context(parameterList);

    auto assemblyOptions = context.assemblyOptions();
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
      KernelFunctor;
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
  shared_ptr<ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>> newOp(new Op(domain, range, dualToRange, label, symmetry,
                              KernelFunctor(), TransformationFunctor(),
                              TransformationFunctor(), integral));
  return newOp;

}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
{

    Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<KernelType>
      KernelFunctor;
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

  shared_ptr<Op> newOp(new Op(domain, range, dualToRange, label, symmetry,
                              KernelFunctor(), TransformationFunctor(),
                              TransformationFunctor(), integral));

  return newOp;
}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceAdjointDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
{
    Context<BasisFunctionType, ResultType> context(parameterList);

  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Laplace3dAdjointDoubleLayerPotentialKernelFunctor<KernelType>
      KernelFunctor;
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
  shared_ptr<Op> newOp(new Op(domain, range, dualToRange, label, symmetry,
                              KernelFunctor(), TransformationFunctor(),
                              TransformationFunctor(), integral));
  return newOp;

}

template<typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceHypersingularBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
{

  Context<BasisFunctionType, ResultType> context(parameterList);
  const AssemblyOptions &assemblyOptions = context.assemblyOptions();

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
      KernelFunctor;
  typedef Fiber::SurfaceCurl3dFunctor<CoordinateType> TransformationFunctor;
  typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctorExt<
      BasisFunctionType, KernelType, ResultType, 3> IntegrandFunctor;

  typedef GeneralElementarySingularIntegralOperator<BasisFunctionType, KernelType,
                                               ResultType> Op;

  shared_ptr<Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                            ResultType>> integral;

  if (shouldUseBlasInQuadrature(assemblyOptions, *domain, *dualToRange)) {
    integral.reset(new Fiber::TypicalTestScalarKernelTrialIntegral<
        BasisFunctionType, KernelType, ResultType>());
  } else {
    integral.reset(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
        IntegrandFunctor()));
  }

  shared_ptr<Op> newOp(
      new Op(domain, range, dualToRange, label, symmetry, KernelFunctor(),
             TransformationFunctor(), TransformationFunctor(), integral));
  return newOp;

}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>> laplaceSingleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        const ParameterList& parameterList)
{

        typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
        shared_ptr<Matrix<CoordinateType>> pointsPtr(
                new Matrix<CoordinateType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, ResultType>> op(
                new Laplace3dSingleLayerPotentialOperator<BasisFunctionType, ResultType>());
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>> laplaceDoubleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        const ParameterList& parameterList)
{

        typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
        shared_ptr<Matrix<CoordinateType>> pointsPtr(
                new Matrix<CoordinateType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, ResultType>> op(
                new Laplace3dDoubleLayerPotentialOperator<BasisFunctionType, ResultType>());
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();
}

}

#endif

