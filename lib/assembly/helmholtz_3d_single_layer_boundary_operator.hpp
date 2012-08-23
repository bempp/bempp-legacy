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

#ifndef bempp_helmholtz_3d_single_layer_boundary_operator_hpp
#define bempp_helmholtz_3d_single_layer_boundary_operator_hpp

#include "helmholtz_3d_boundary_operator_base.hpp"
#include "boundary_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
struct Helmholtz3dSingleLayerBoundaryOperatorImpl;

/** \ingroup helmholtz_3d
 *  \brief Single-layer-potential boundary operator for the Helmholtz equation in 3D.
 *
 *  \tparam BasisFunctionType
 *    Type of the values of the basis functions into which functions acted upon
 *    by the operator are expanded. It can take the following values: \c float,
 *    \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename BasisFunctionType_>
class Helmholtz3dSingleLayerBoundaryOperator :
        public Helmholtz3dBoundaryOperatorBase<
                Helmholtz3dSingleLayerBoundaryOperatorImpl<BasisFunctionType_>,
                BasisFunctionType_>
{
    typedef Helmholtz3dBoundaryOperatorBase<
    Helmholtz3dSingleLayerBoundaryOperatorImpl<BasisFunctionType_>,
    BasisFunctionType_> Base;
public:
    /** \copydoc ElementaryIntegralOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc ElementaryIntegralOperator::KernelType */
    typedef typename Base::KernelType KernelType;
    /** \copydoc ElementaryIntegralOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryIntegralOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryIntegralOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryIntegralOperator::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ElementaryIntegralOperator::TestKernelTrialIntegral */
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    /** \copydoc Helmholtz3dBoundaryOperatorBase::Helmholtz3dBoundaryOperatorBase */
    Helmholtz3dSingleLayerBoundaryOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            KernelType waveNumber,
            const std::string& label = "");
};

/** \relates Helmholtz3dSingleLayerBoundaryOperator
 *  \brief Construct a BoundaryOperator object wrapping a
 *  Helmholtz3dSingleLayerBoundaryOperator.
 *
 *  This is a convenience function that creates a
 *  Helmholtz3dSingleLayerBoundaryOperator, immediately wraps it in a
 *  BoundaryOperator and returns the latter object.
 *
 *  \param[in] context
 *    A Context object that will be used to build the weak form of the
 *    boundary operator when necessary.
 *  \param[in] domain
 *    Function space being the domain of the boundary operator.
 *  \param[in] range
 *    Function space being the range of the boundary operator.
 *  \param[in] dualToRange
 *    Function space dual to the the range of the boundary operator.
 *  \param[in] waveNumber
 *    Wave number. See \ref helmholtz_3d for its definition.
 *  \param[in] label
 *    Textual label of the operator (optional, used for debugging).
 *
 *  None of the shared pointers may be null and the spaces \p range and \p
 *  dualToRange must be defined on the same grid, otherwise an exception is
 *  thrown. */
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename Helmholtz3dSingleLayerBoundaryOperator<BasisFunctionType>::ResultType>
helmholtz3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename Helmholtz3dSingleLayerBoundaryOperator<BasisFunctionType>::ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename Helmholtz3dSingleLayerBoundaryOperator<BasisFunctionType>::KernelType waveNumber,
        const std::string& label = "");

} // namespace Bempp

#endif
