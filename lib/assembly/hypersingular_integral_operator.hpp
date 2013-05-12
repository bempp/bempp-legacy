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
/** \endcond */

/** \ingroup abstract_boundary_operators
 *  \brief Hypersingular integral boundary operator.
 *
 *  TODO: write documentation.
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
    typedef Fiber::LocalAssemblerForOperators<ResultType> LocalAssembler;
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

    /** \brief Return the collection of kernel functions occurring in the
     *  weak form of this operator. */
    virtual const CollectionOfKernels& offDiagonalKernels() const = 0;

    /** \brief Return the collection of test-function transformations occurring
     *  in the weak form of this operator. */
    virtual const CollectionOfBasisTransformations&
    offDiagonalTestTransformations() const = 0;

    /** \brief Return the collection of trial-function transformations occurring
     *  in the weak form of this operator. */
    virtual const CollectionOfBasisTransformations&
    offDiagonalTrialTransformations() const = 0;

    /** \brief Return an object representing the integral that is the weak form
     *  of this operator.
     *
     *  Subclasses of #TestKernelTrialIntegral implement functions that evaluate
     *  the integral using the data provided by a #CollectionOfKernels
     *  representing the kernel functions occurring in the integrand and a pair
     *  of #CollectionOfBasisTransformations objects representing the test and
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
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
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
