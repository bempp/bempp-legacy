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

#ifndef bempp_modified_helmholtz_3d_double_layer_potential_operator_hpp
#define bempp_modified_helmholtz_3d_double_layer_potential_operator_hpp

#include "modified_helmholtz_3d_potential_operator_base.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct ModifiedHelmholtz3dDoubleLayerPotentialOperatorImpl;
/** \endcond */

/** \ingroup helmholtz_3d
 *  \brief Double-layer potential operator for the modified Helmholtz equation in 3D.
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into which functions acted upon
 *    by the operator are expanded. It can take the following values: \c float,
 *    \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename BasisFunctionType_>
class ModifiedHelmholtz3dDoubleLayerPotentialOperator :
        public ModifiedHelmholtz3dPotentialOperatorBase<
        ModifiedHelmholtz3dDoubleLayerPotentialOperatorImpl<BasisFunctionType_>,
        BasisFunctionType_>
{
    typedef ModifiedHelmholtz3dPotentialOperatorBase<
    ModifiedHelmholtz3dDoubleLayerPotentialOperatorImpl<BasisFunctionType_>,
    BasisFunctionType_>
    Base;
public:
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::KernelType */
    typedef typename Base::KernelType KernelType;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::KernelTrialIntegral */
    typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::ModifiedHelmholtz3dPotentialOperatorBase */
    ModifiedHelmholtz3dDoubleLayerPotentialOperator(KernelType waveNumber);
    /** \copydoc ModifiedHelmholtz3dPotentialOperatorBase::~ModifiedHelmholtz3dPotentialOperatorBase */
    virtual ~ModifiedHelmholtz3dDoubleLayerPotentialOperator();
};

} // namespace Bempp

#endif
