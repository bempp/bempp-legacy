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

#ifndef bempp_elementary_local_operator_hpp
#define bempp_elementary_local_operator_hpp

#include "../common/common.hpp"

#include "elementary_abstract_boundary_operator.hpp"

#include "abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
template <typename BasisFunctionType, typename ResultType> class TestTrialIntegral;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
/** \endcond */

/** \ingroup identity
 *  \brief Abstract base class of local elementary operators.
 *
 *  See AbstractBoundaryOperator for the documentation of the template
 *  parameters.
 */
template <typename BasisFunctionType_, typename ResultType_>
class ElementaryLocalOperator :
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
    /** \brief Type of the appropriate instantiation of
     *  Fiber::CollectionOfBasisTransformations. */
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    /** \brief Type of the appropriate instantiation of Fiber::TestTrialIntegral. */
    typedef Fiber::TestTrialIntegral<BasisFunctionType, ResultType>
    TestTrialIntegral;

    /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
    ElementaryLocalOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label,
            int symmetry);

    /** \brief Return true. */
    virtual bool isLocal() const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:    
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
     *  Subclasses of #TestTrialIntegral implement functions that evaluate
     *  the integral using the data provided by a pair
     *  of #CollectionOfBasisTransformations objects representing the test and
     *  trial basis function transformations occurring in the integrand. */
    virtual const TestTrialIntegral& integral() const = 0;

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

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInSparseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

} // namespace Bempp

#endif
