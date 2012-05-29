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

#ifndef bempp_modified_helmholtz_3d_single_layer_potential_hpp
#define bempp_modified_helmholtz_3d_single_layer_potential_hpp

#include "elementary_weakly_singular_integral_operator.hpp"
#include "../common/scalar_traits.hpp"
#include "../fiber/modified_helmholtz_3d_single_layer_potential_kernel.hpp"
#include "../fiber/scalar_function_value.hpp"

namespace Bempp
{

/** \ingroup modified_helmholtz_3d
 *  \brief Adjoint double-layer-potential operator for the modified Helmholtz
 *  equation in 3D.
 *
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions.
 *  \tparam KernelType
 *    Type used to represent the values of the kernel.
 *  \tparam ResultType
 *    Type used to represent entries in the discrete form of the operator.
 *
 *  All three template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>. All
 *  types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType
 *  is by default set to "larger" of \p BasisFunctionType and \p KernelType,
 *  e.g. for \p BasisFunctionType = \c double and \p KernelType =
 *  <tt>std::complex<double></tt> it is set to <tt>std::complex<double></tt>.
 *  You should override that only if you set both \p BasisFunctionType and \p
 *  KernelType to a real type, but you want the entries of the operator's weak
 *  form to be stored as complex numbers.
 *
 *  Note that setting \p KernelType to a real type implies that the wave number
 *  must also be chosen purely real.
 *
 *  \see modified_helmholtz_3d
 */
template <typename BasisFunctionType, typename KernelType,
          typename ResultType = typename Coercion<BasisFunctionType, KernelType>::Type>
class ModifiedHelmholtz3dSingleLayerPotential :
        public ElementaryWeaklySingularIntegralOperator<
        BasisFunctionType, KernelType, ResultType>
{
    typedef ElementaryWeaklySingularIntegralOperator<
    BasisFunctionType, KernelType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct the operator.
     *
     * \param testSpace Test function space.
     * \param trialSpace Trial function space.
     * \param waveNumber Wave number.
     *
     * See \ref modified_helmholtz_3d for the definition of the wave number. */
    ModifiedHelmholtz3dSingleLayerPotential(
            const Space<BasisFunctionType>& testSpace,
            const Space<BasisFunctionType>& trialSpace,
            KernelType waveNumber);

private:
    virtual const Fiber::Kernel<KernelType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::Expression<CoordinateType>& testExpression() const {
        return m_expression;
    }

    virtual const Fiber::Expression<CoordinateType>& trialExpression() const {
        return m_expression;
    }

private:
    Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernel<KernelType> m_kernel;
    Fiber::ScalarFunctionValue<CoordinateType> m_expression;
};

} // namespace Bempp

#endif
