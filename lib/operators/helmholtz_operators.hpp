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

#ifndef bempp_helmholtz_operators_hpp
#define bempp_helmholtz_operators_hpp

#include "../assembly/symmetry.hpp"
#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../operators/modified_helmholtz_operators.hpp"
#include "helmholtz_3d_far_field_double_layer_potential_operator.hpp"
#include "helmholtz_3d_far_field_single_layer_potential_operator.hpp"
#include "../fiber/scalar_traits.hpp"

#include <complex>

namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType> class ElementaryIntegralOperator;
template <typename BasisFunctionType, typename ResultType> class PotentialOperator;
template <typename BasisFunctionType> class Space;

template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, 
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
helmholtzSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY)
{

    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;


    return modifiedHelmholtzSingleLayerBoundaryOperator<BasisFunctionType, KernelType, ResultType>(
            parameterList, domain, range, dualToRange, waveNumber/KernelType(0,1), label, symmetry);


}

template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, 
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
helmholtzDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY)
{

    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;


    return modifiedHelmholtzDoubleLayerBoundaryOperator<BasisFunctionType, KernelType, ResultType>(
            parameterList, domain, range, dualToRange, waveNumber/KernelType(0,1), label, symmetry);


}

template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, 
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
helmholtzAdjointDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY)
{

    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;


    return modifiedHelmholtzAdjointDoubleLayerBoundaryOperator<BasisFunctionType, KernelType, ResultType>(
            parameterList, domain, range, dualToRange, waveNumber/KernelType(0,1), label, symmetry);


}

template<typename BasisFunctionType>
shared_ptr<const ElementaryIntegralOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType, 
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>>
helmholtzHypersingularBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY)
{

    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;


    return modifiedHelmholtzHypersingularBoundaryOperator<BasisFunctionType, KernelType, ResultType>(
            parameterList, domain, range, dualToRange, waveNumber/KernelType(0,1), label, symmetry);


}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> helmholtzSingleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;

    return modifiedHelmholtzSingleLayerPotentialOperator(space, evaluationPoints, waveNumber/KernelType(0,1), parameterList);

}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> helmholtzDoubleLayerPotentialOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType ResultType;

    return modifiedHelmholtzDoubleLayerPotentialOperator(space, evaluationPoints, waveNumber/KernelType(0,1), parameterList);

}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> helmholtzSingleLayerFarFieldOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<BasisFunctionType>> pointsPtr(
                new Matrix<BasisFunctionType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> op(
                new Helmholtz3dFarFieldSingleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();

}

template <typename BasisFunctionType>
shared_ptr<const DiscreteBoundaryOperator<typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> helmholtzDoubleLayerFarFieldOperator(
        const shared_ptr<const Space<BasisFunctionType>>& space,
        const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>& evaluationPoints,
        typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const ParameterList& parameterList)
{

        shared_ptr<Matrix<BasisFunctionType>> pointsPtr(
                new Matrix<BasisFunctionType>(evaluationPoints));

        shared_ptr<PotentialOperator<BasisFunctionType, typename Fiber::ScalarTraits<BasisFunctionType>::ComplexType>> op(
                new Helmholtz3dFarFieldDoubleLayerPotentialOperator<double>(waveNumber));
        return op->assemble(space, pointsPtr, parameterList).discreteOperator();

}

}

#endif
