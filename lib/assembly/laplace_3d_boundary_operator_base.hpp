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

#include "abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Bempp
{

template <typename Impl, typename BasisFunctionType, typename ResultType>
class Laplace3dBoundaryOperatorBase;

template <typename BasisFunctionType>
class Laplace3dBoundaryOperatorId : public AbstractBoundaryOperatorId
{
public:
    template <typename Impl, typename ResultType>
    explicit Laplace3dBoundaryOperatorId(
            const Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual void dump() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const std::type_info& m_typeInfo;
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

/** \ingroup laplace_3d
 *  \brief Base class for boundary operators for the Laplace equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation pointer.
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions.
 *  \tparam ResultType
 *    Type used to represent entries in the discrete form of the operator.
 *
 *  The latter two template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType
 *  is by default set to \p BasisFunctionType. You should override that only if
 *  you set \p BasisFunctionType to a real type, but you want the entries of
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
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::KernelType KernelType;
    typedef typename Base::ResultType ResultType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    Laplace3dBoundaryOperatorBase(
            const Space<BasisFunctionType>& domain,
            const Space<BasisFunctionType>& range,
            const Space<BasisFunctionType>& dualToRange,
            const std::string& label = "");
    Laplace3dBoundaryOperatorBase(
            const Laplace3dBoundaryOperatorBase& other);
    virtual ~Laplace3dBoundaryOperatorBase();

    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

private:
    virtual const CollectionOfKernels& kernels() const;
    virtual const CollectionOfBasisTransformations&
    testTransformations() const;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;
    virtual const TestKernelTrialIntegral& integral() const;

private:
    boost::scoped_ptr<Impl> m_impl;
    boost::shared_ptr<AbstractBoundaryOperatorId> m_id;
};

} // namespace Bempp

#endif
