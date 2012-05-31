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

#ifndef bempp_helmholtz_3d_hypersingular_operator_hpp
#define bempp_helmholtz_3d_hypersingular_operator_hpp

#include "elementary_singular_integral_operator.hpp"
#include "../common/scalar_traits.hpp"
#include "../fiber/expression_list.hpp"
#include "../fiber/modified_helmholtz_3d_single_layer_potential_kernel.hpp"
#include "../fiber/scalar_function_value_times_normal_3d.hpp"
#include "../fiber/surface_curl_3d.hpp"

namespace Bempp
{

/** \ingroup helmholtz_3d
 *  \brief Hypersingular operator for the Helmholtz equation in 3D.
 *
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions. It can take the
 *    following values: \c float, \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename BasisFunctionType>
class Helmholtz3dHypersingularOperator :
        public ElementarySingularIntegralOperator<
        BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType,
        typename ScalarTraits<BasisFunctionType>::ComplexType>
{
public:
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef KernelType ResultType;
private:
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType, KernelType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct the operator.
     *
     * \param testSpace Test function space.
     * \param trialSpace Trial function space.
     * \param waveNumber Wave number.
     *
     * See \ref helmholtz_3d for the definition of the wave number. */
    Helmholtz3dHypersingularOperator(const Space<BasisFunctionType>& testSpace,
                                     const Space<BasisFunctionType>& trialSpace,
                                     KernelType waveNumber);

private:
    virtual const Fiber::Kernel<KernelType>& kernel() const {
        return m_kernel;
    }

    virtual const Fiber::ExpressionList<ResultType>& testExpressionList() const {
        return m_expressionList;
    }

    virtual const Fiber::ExpressionList<ResultType>& trialExpressionList() const {
        return m_expressionList;
    }

private:
    Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernel<KernelType> m_kernel;
    Fiber::SurfaceCurl3d<CoordinateType> m_surfaceCurl;
    Fiber::ScalarFunctionValueTimesNormal3d<CoordinateType> m_valueTimesNormal;
    Fiber::ExpressionList<ResultType> m_expressionList;
};

} // namespace Bempp

#endif
