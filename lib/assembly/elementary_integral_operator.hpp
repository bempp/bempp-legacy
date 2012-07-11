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

template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename KernelType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;
template <typename ResultType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

class EvaluationOptions;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ValueType> class InterpolatedFunction;
template <typename BasisFunctionType, typename ResultType> class WeakFormAcaAssemblyHelper;

/** \ingroup assembly
 *  \brief Elementary integral operator.
 */
template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class ElementaryIntegralOperator :
        public ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::ResultType ResultType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    typedef typename Base::LocalAssembler LocalAssembler;
    typedef KernelType_ KernelType;
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
    typedef Fiber::TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
    TestKernelTrialIntegral;

    ElementaryIntegralOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                               const shared_ptr<const Space<BasisFunctionType> >& range,
                               const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                               const std::string& label = "",
                               const Symmetry = NO_SYMMETRY);

    virtual int trialComponentCount() const;
    virtual int testComponentCount() const;

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

    virtual bool isRegular() const = 0;

protected:
    /** @}
        \name Weak form assembly
        @{ */
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    virtual const CollectionOfKernels& kernels() const = 0;
    virtual const CollectionOfBasisTransformations&
    testTransformations() const = 0;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const = 0;
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
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const;

    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInAcaMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
    /** @} */

};

} // namespace Bempp

#endif
