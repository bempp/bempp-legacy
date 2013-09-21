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

#ifndef bempp_hypersingular_integral_operator_hpp
#define bempp_hypersingular_integral_operator_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator.hpp"
#include "../common/types.hpp"
#include "../fiber/types.hpp"

#include <vector>
#include "../common/armadillo_fwd.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfShapesetTransformations;
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
/** \endcond */

/** \ingroup abstract_boundary_operators
 *  \brief Hypersingular integral boundary operator.
 *
 *  This abstract class represents an integral boundary operator \f$\mathcal
 *  A\f$ providing two alternative ways to evaluate its weak form \f$\langle
 *  \phi, \mathcal A \psi \rangle\f$. One of them must be valid for any pair of
 *  test and trial functions, \f$\phi\f$ and \f$\psi\f$, and it can be
 *  expressed by a formula of an arbitrary form; in particular, it can contain
 *  differential operators, e.g. surface div and curl, acting on the test and
 *  trial functions. The other one needs only be valid for test and trial
 *  functions with nonoverlapping support, and it must be of the form
 *
 *  \f[
 *      \int_\Gamma \int_\Sigma
 *      \phi^*(x) \, K(x, y) \, \psi(y) \,
 *      \mathrm{d}\Gamma(x) \,\mathrm{d}\Sigma(y),
 *  \f]
 *
 *  where \f$K(x, y)\f$ is an arbitrary kernel; note that the formula contains
 *  "bare" test and trial functions, without any operators acting on them.
 *
 *  When the discrete weak form of an integral operator represented by a
 *  subclass of this class is assembled in the hybrid ACA mode (\see
 *  AcaOptions::mode), the first representation of the weak form is used in the
 *  assembly of nonadmissible blocks and the second in the assembly of
 *  admissible ("off-diagonal") blocks.
 *
 *  The first representation of the weak form is determined by the
 *  implementation of the virtual functions integral(), kernels(),
 *  testTransformations() and trialTransformations(). The second representation
 *  is determined by the implementation of the virtual functions
 *  offDiagonalIntegral(), offDiagonalKernels(),
 *  offDiagonalTestTransformations() and offDiagonalTrialTransformations().
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
class HypersingularIntegralOperator :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    /** \brief Type of the appropriate instantiation of Fiber::LocalAssemblerForOperators. */
    typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
    /** \brief Type of the values of the (components of the) kernel functions. */
    typedef KernelType_ KernelType;
    /** \brief Type of the appropriate instantiation of Fiber::CollectionOfShapesetTransformations. */
    typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
    CollectionOfShapesetTransformations;
    /** \brief Type of the appropriate instantiation of Fiber::CollectionOfBasisTransformations.
     *
     *  \deprecated This type is deprecated; use CollectionOfShapesetTransformations
     *  instead. */
    typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    /** \brief Type of the appropriate instantiation of Fiber::CollectionOfKernels. */
    typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
    /** \brief Type of the appropriate instantiation of Fiber::TestKernelTrialIntegral. */
    typedef Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
    TestKernelTrialIntegral;

    /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
    HypersingularIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label,
            int symmetry);

    /** \brief Return false. */
    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    /** \brief Return the collection of kernel functions occurring in the
     *  weak form of this operator. */
    virtual const CollectionOfKernels& kernels() const = 0;

    /** \brief Return the collection of test-function transformations occurring
     *  in the first representation of the weak form of this operator (valid for
     *  all pairs of test and trial functions). */
    virtual const CollectionOfShapesetTransformations&
    testTransformations() const = 0;

    /** \brief Return the collection of trial-function transformations occurring
     *  in the first representation of the weak form of this operator (valid for
     *  all pairs of test and trial functions). */
    virtual const CollectionOfShapesetTransformations&
    trialTransformations() const = 0;

    /** \brief Return an object representing the integral that
     *  is the first representation of the weak form of this operator (valid for
     *  all pairs of test and trial functions).
     *
     *  Subclasses of #TestKernelTrialIntegral implement functions that evaluate
     *  the integral using the data provided by a #CollectionOfKernels
     *  representing the kernel functions occurring in the integrand and a pair
     *  of #CollectionOfShapesetTransformations objects representing the test and
     *  trial function transformations occurring in the integrand. */
    virtual const TestKernelTrialIntegral& integral() const = 0;

    /** \brief Return the collection of kernel functions occurring in the
     *  in the second representation of the weak form of this operator (valid at
     *  least for pairs of test and trial functions with nonoverlapping
     *  supports). */
    virtual const CollectionOfKernels& offDiagonalKernels() const = 0;

    /** \brief Return the collection of test-function transformations occurring
     *  in the second representation of the weak form of this operator (valid at
     *  least for pairs of test and trial functions with nonoverlapping
     *  supports).
     *
     *  It should normally be the mapping of shape functions (defined on the
     *  reference element) to basis functions (defined on the physical
     *  elements). */
    virtual const CollectionOfShapesetTransformations&
    offDiagonalTestTransformations() const = 0;

    /** \brief Return the collection of trial-function transformations occurring
     *  in the second representation of the weak form of this operator (valid at
     *  least for pairs of test and trial functions with nonoverlapping
     *  supports).
     *
     *  It should normally be the mapping of shape functions (defined on the
     *  reference element) to basis functions (defined on the physical
     *  elements). */
    virtual const CollectionOfShapesetTransformations&
    offDiagonalTrialTransformations() const = 0;

    /** \brief Return an object representing the integral that
     *  is the second representation of the weak form of this operator (valid at
     *  least for pairs of test and trial functions with nonoverlapping
     *  supports).
     *
     *  Subclasses of #TestKernelTrialIntegral implement functions that evaluate
     *  the integral using the data provided by a #CollectionOfKernels
     *  representing the kernel functions occurring in the integrand and a pair
     *  of #CollectionOfShapesetTransformations objects representing the test and
     *  trial basis function transformations occurring in the integrand. */
    virtual const TestKernelTrialIntegral& offDiagonalIntegral() const = 0;

private:
    /** \cond PRIVATE */

    std::pair<shared_ptr<LocalAssembler>, shared_ptr<LocalAssembler> >
    makeAssemblers(
            const QuadratureStrategy& quadStrategy,
            const AssemblyOptions& options) const;

    std::pair<shared_ptr<LocalAssembler>, shared_ptr<LocalAssembler> >
    reallyMakeAssemblers(
            const QuadratureStrategy& quadStrategy,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Shapeset<BasisFunctionType>*> >& testShapesets,
            const shared_ptr<const std::vector<const Fiber::Shapeset<BasisFunctionType>*> >& trialShapesets,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            VerbosityLevel::Level verbosityLevel,
            bool cacheSingularIntegrals,
            bool makeSeparateOffDiagonalAssembler) const;

    shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInternal(
            LocalAssembler& standardAssembler, LocalAssembler& offDiagonalAssembler,
            const Context<BasisFunctionType, ResultType>& context) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const Context<BasisFunctionType, ResultType>& context) const;
    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInAcaMode(
            LocalAssembler& standardAssembler, LocalAssembler& offDiagonalAssembler,
            const Context<BasisFunctionType, ResultType>& context) const;

    /** \endcond */
};

} // namespace Bempp

#endif
