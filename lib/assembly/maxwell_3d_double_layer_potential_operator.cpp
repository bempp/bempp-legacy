// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "maxwell_3d_double_layer_potential_operator.hpp"
#include "helmholtz_3d_potential_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_maxwell_3d_double_layer_operators_kernel_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_potential_operator_integrand_functor.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_kernel_trial_integral.hpp"

namespace Bempp {

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct Maxwell3dDoubleLayerPotentialOperatorImpl {
  typedef Maxwell3dDoubleLayerPotentialOperatorImpl<BasisFunctionType> This;
  typedef Helmholtz3dPotentialOperatorBase<This, BasisFunctionType>
      PotentialOperatorBase;
  typedef typename PotentialOperatorBase::KernelType KernelType;
  typedef typename PotentialOperatorBase::ResultType ResultType;
  typedef typename PotentialOperatorBase::CoordinateType CoordinateType;

  typedef Fiber::ModifiedMaxwell3dDoubleLayerOperatorsKernelFunctor<KernelType>
      KernelFunctor;
  typedef Fiber::HdivFunctionValueFunctor<CoordinateType> TransformationFunctor;
  typedef Fiber::ModifiedMaxwell3dDoubleLayerPotentialOperatorIntegrandFunctor<
      BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

  Maxwell3dDoubleLayerPotentialOperatorImpl(KernelType waveNumber)
      : kernels(KernelFunctor(waveNumber / KernelType(0., 1.))),
        transformations(TransformationFunctor()), integral(IntegrandFunctor()) {
  }

  Fiber::DefaultCollectionOfKernels<KernelFunctor> kernels;
  Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
      transformations;
  Fiber::DefaultKernelTrialIntegral<IntegrandFunctor> integral;
};
/** \endcond */

template <typename BasisFunctionType>
Maxwell3dDoubleLayerPotentialOperator<BasisFunctionType>::
    Maxwell3dDoubleLayerPotentialOperator(KernelType waveNumber)
    : Base(waveNumber) {}

template <typename BasisFunctionType>
Maxwell3dDoubleLayerPotentialOperator<
    BasisFunctionType>::~Maxwell3dDoubleLayerPotentialOperator() {}

#define INSTANTIATE_BASE_HELMHOLTZ_DOUBLE_POTENTIAL(BASIS)                     \
  template class Helmholtz3dPotentialOperatorBase<                             \
      Maxwell3dDoubleLayerPotentialOperatorImpl<BASIS>, BASIS>
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_BASE_HELMHOLTZ_DOUBLE_POTENTIAL);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    Maxwell3dDoubleLayerPotentialOperator);

} // namespace Bempp
