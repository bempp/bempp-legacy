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

#ifndef bempp_laplace_3d_potential_operator_base_hpp
#define bempp_laplace_3d_potential_operator_base_hpp

#include "elementary_potential_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup laplace_3d
 *  \brief Base class for potential operators related to the Laplace equation in
 *  3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation object.
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType_
 *    Type of the values of the potential.
 *
 *  The latter two template parameters can take the following values: \c float,
 *  \c double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType_
 *  is by default set to \p BasisFunctionType_. You should override that only if
 *  you set \p BasisFunctionType_ to a real type, but you want the potential
 *  values to be stored as complex numbers.
 *
 *  \see laplace_3d */
template <typename Impl,
          typename BasisFunctionType_, typename ResultType_ = BasisFunctionType_>
class Laplace3dPotentialOperatorBase :
        public ElementaryPotentialOperator<
        BasisFunctionType_,
        typename ScalarTraits<ResultType_>::RealType,
        ResultType_>
{
    typedef ElementaryPotentialOperator<
    BasisFunctionType_,
    typename ScalarTraits<ResultType_>::RealType,
    ResultType_>
    Base;
public:
    /** \brief Type of the values of the basis functions into
     *  which functions acted upon by the operator are expanded. */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \brief Type of the values of kernel functions. */
    typedef typename Base::KernelType KernelType;
    /** \brief Type of the values of the potential. */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryPotentialOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryPotentialOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryPotentialOperator::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ElementaryPotentialOperator::KernelTrialIntegral */
    typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

    /** \brief Constructor. */
    Laplace3dPotentialOperatorBase();
    /** \brief Copy constructor. */
    Laplace3dPotentialOperatorBase(const Laplace3dPotentialOperatorBase& other);
    /** \brief Destructor. */
    virtual ~Laplace3dPotentialOperatorBase();

private:
    virtual const CollectionOfKernels& kernels() const;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;
    virtual const KernelTrialIntegral& integral() const;

private:
    /** \cond PRIVATE */
    boost::scoped_ptr<Impl> m_impl;
    /** \endcond */
};

} // namespace Bempp

#endif
