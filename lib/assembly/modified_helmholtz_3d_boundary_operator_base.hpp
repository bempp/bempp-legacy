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

#ifndef bempp_modified_helmholtz_3d_boundary_operator_base_hpp
#define bempp_modified_helmholtz_3d_boundary_operator_base_hpp

#include "elementary_singular_integral_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup helmholtz_3d
 *  \brief Base class for boundary operators for the modified Helmholtz
 *  equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation pointer.
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions.
 *  \tparam KernelType
 *    Type used to represent the values of the kernel.
 *  \tparam ResultType
 *    Type used to represent entries in the discrete form of the operator.
 *
 *  The latter three template parameters can take the following values: \c
 *  float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. All types must have the same precision: for
 *  instance, mixing \c float with <tt>std::complex<double></tt> is not
 *  allowed. The parameter \p ResultType is by default set to "larger" of \p
 *  BasisFunctionType and \p KernelType, e.g. for \p BasisFunctionType = \c
 *  double and \p KernelType = <tt>std::complex<double></tt> it is set to
 *  <tt>std::complex<double></tt>. You should override that only if you set
 *  both \p BasisFunctionType and \p KernelType to a real type, but you want
 *  the entries of the operator's weak form to be stored as complex numbers.
 *
 *  Note that setting \p KernelType to a real type implies that the wave number
 *  must also be chosen purely real.
 *
 *  \see modified_helmholtz_3d
 */
template <typename Impl, typename BasisFunctionType_,
          typename KernelType_, typename ResultType_>
class ModifiedHelmholtz3dBoundaryOperatorBase :
        public ElementarySingularIntegralOperator<
        BasisFunctionType_, KernelType_, ResultType_>
{
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType_, KernelType_, ResultType_>
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

    ModifiedHelmholtz3dBoundaryOperatorBase(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            KernelType waveNumber,
            const std::string& label = "");
    ModifiedHelmholtz3dBoundaryOperatorBase(
            const ModifiedHelmholtz3dBoundaryOperatorBase& other);
    virtual ~ModifiedHelmholtz3dBoundaryOperatorBase();

    void setWaveNumber(KernelType waveNumber);
    KernelType waveNumber() const;

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
    shared_ptr<AbstractBoundaryOperatorId> m_id;
};

} // namespace Bempp

#endif
