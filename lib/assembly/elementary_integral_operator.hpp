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

#ifndef bempp_elementary_integral_operator_hpp
#define bempp_elementary_integral_operator_hpp

#include "../common/common.hpp"

#include "elementary_abstract_boundary_operator.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/types.hpp"
#include "../fiber/types.hpp"

#include <vector>
#include "../common/armadillo_fwd.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename KernelType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class EvaluationOptions;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ValueType> class InterpolatedFunction;
template <typename BasisFunctionType, typename ResultType> class WeakFormAcaAssemblyHelper;
/** \endcond */

/** \ingroup abstract_boundary_operators
 *  \brief Elementary integral boundary operator.
 *
 *  This class represents an integral boundary operator \f$\mathcal A\f$ whose
 *  weak form is
 *  \f[
 *    \langle \phi, \mathcal A \psi \rangle \equiv
 *    \int_\Gamma \int_\Gamma F[\phi(x), \psi(y)] \, d\Gamma(x) \, d\Gamma(y),
 *  \f]
 *  where \f$\Gamma\f$ is a surface embedded in a space of dimension higher by
 *  one and the integrand \f$F[\phi(x), \psi(y)]\f$ is an arbitrary bilinear (or
 *  sesquilinear) form of the *test function* \f$\phi\f$ belonging to the space
 *  dual to the range of \f$\mathcal A\f$ and the *trial function* \f$\psi\f$
 *  belonging to the domain of \f$\mathcal A\f$. In the simplest and most
 *  common case, \f$F[\phi(x), \psi(y)]\f$ is just
 *  \f[
 *    F[\phi(x), \psi(y)] = \phi^*(x) K(x, y) \psi(y),
 *  \f]
 *  where \f$K(x, y)\f$ is a *kernel function* and the asterisk denotes complex
 *  conjugation. For more complex operators, \f$F\f$ might involve some
 *  transformations of the test and trial functions (e.g. their surface
 *  divergence or curl), the kernel function might be a tensor, or \f$F\f$
 *  might consist of several terms. The exact form of \f$F\f$ for a particular
 *  boundary operator is determined by the implementation of the virtual
 *  functions integral(), kernels(), testTransformations() and
 *  trialTransformations().
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam KernelType_
 *    Type of the values of the (components of the) kernel functions occurring
 *    in the integrand of the operator.
 *  \tparam ResultType_
 *    Type used to represent elements of the weak form of the operator.
 *
 *  All three template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>. All
 *  types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If either \p
 *  BasisFunctionType_ or \p KernelType_ is a complex type, then \p ResultType_
 *  must be set to the same type. */
template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class ElementaryIntegralOperator :
        public ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc ElementaryAbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc ElementaryAbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryAbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryAbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    /** \copydoc ElementaryAbstractBoundaryOperator::LocalAssembler */
    typedef typename Base::LocalAssembler LocalAssembler;
    /** \brief Type of the values of the (components of the) kernel functions. */
    typedef KernelType_ KernelType;
    /** \brief Type of the appropriate instantiation of Fiber::CollectionOfBasisTransformations. */
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    /** \brief Type of the appropriate instantiation of Fiber::CollectionOfKernels. */
    typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
    /** \brief Type of the appropriate instantiation of Fiber::TestKernelTrialIntegral. */
    typedef Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
    TestKernelTrialIntegral;

    /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
    ElementaryIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label,
            int symmetry);

    /** \brief Return false. */
    virtual bool isLocal() const;

    /** \brief Return whether applying this operator to a regular function
     *  yields a regular integral. */
    virtual bool isRegular() const = 0;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    /** \brief Return the collection of kernel functions occurring in the
     *  weak form of this operator. */
    virtual const CollectionOfKernels& kernels() const = 0;

    /** \brief Return the collection of test-function transformations occurring
     *  in the weak form of this operator. */
    virtual const CollectionOfBasisTransformations&
    testTransformations() const = 0;

    /** \brief Return the collection of trial-function transformations occurring
     *  in the weak form of this operator. */
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const = 0;

    /** \brief Return an object representing the integral that is the weak form
     *  of this operator.
     *
     *  Subclasses of #TestKernelTrialIntegral implement functions that evaluate
     *  the integral using the data provided by a #CollectionOfKernels
     *  representing the kernel functions occurring in the integrand and a pair
     *  of #CollectionOfBasisTransformations objects representing the test and
     *  trial basis function transformations occurring in the integrand. */
    virtual const TestKernelTrialIntegral& integral() const = 0;

    virtual std::auto_ptr<LocalAssembler> makeAssemblerImpl(
            const QuadratureStrategy& quadStrategy,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            VerbosityLevel::Level verbosityLevel,
            bool cacheSingularIntegrals) const;

    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    /** \cond PRIVATE */

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInAcaMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    /** \endcond */
};

} // namespace Bempp

#endif
