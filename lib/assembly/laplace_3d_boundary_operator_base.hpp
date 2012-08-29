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

#ifndef bempp_laplace_3d_boundary_operator_base_hpp
#define bempp_laplace_3d_boundary_operator_base_hpp

#include "elementary_singular_integral_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup laplace_3d
 *  \brief Base class for boundary operators for the Laplace equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation object.
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType_
 *    Type used to represent elements of the weak form form of the operator.
 *
 *  The latter two template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType_
 *  is by default set to \p BasisFunctionType_. You should override that only if
 *  you set \p BasisFunctionType_ to a real type, but you want the entries of
 *  the operator's weak form to be stored as complex numbers.
 *
 *  \see laplace_3d */
template <typename Impl,
          typename BasisFunctionType_, typename ResultType_ = BasisFunctionType_>
class Laplace3dBoundaryOperatorBase :
        public ElementarySingularIntegralOperator<
        BasisFunctionType_,
        typename ScalarTraits<ResultType_>::RealType,
        ResultType_>
{
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType_,
    typename ScalarTraits<ResultType_>::RealType,
    ResultType_>
    Base;
public:
    /** \brief Type of the values of the basis functions into which functions
     *  acted upon by the operator are expanded. */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \brief Type of the values of kernel functions. */
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

    /** \brief Constructor.
     *
     *  \param[in] domain
     *    Function space being the domain of the operator.
     *  \param[in] range
     *    Function space being the range of the operator.
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *  \param[in] label
     *    Textual label of the operator (optional, used for debugging).
     *
     *  None of the shared pointers may be null and the spaces \p range and \p
     *  dualToRange must be defined on the same grid, otherwise an exception is
     *  thrown. */
    Laplace3dBoundaryOperatorBase(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label = "",
            Symmetry symmetry = NO_SYMMETRY);
    Laplace3dBoundaryOperatorBase(
            const Laplace3dBoundaryOperatorBase& other);
    virtual ~Laplace3dBoundaryOperatorBase();

    /** \brief Return the identifier of this operator.
     *
     *  Two boundary operators related to the Laplace equation are treated as
     *  identical, and hence having the same weak form, if they have the same
     *  C++ type (e.g.
     *  <tt>Laplace3dDoubleLayerBoundaryOperator<double, double></tt>),
     *  domain space, range space and space dual to range. */
    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

private:
    virtual const CollectionOfKernels& kernels() const;
    virtual const CollectionOfBasisTransformations&
    testTransformations() const;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;
    virtual const TestKernelTrialIntegral& integral() const;

private:
    /** \cond PRIVATE */
    boost::scoped_ptr<Impl> m_impl;
    shared_ptr<AbstractBoundaryOperatorId> m_id;
    /** \endcond */
};

} // namespace Bempp

#endif
