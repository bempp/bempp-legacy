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


#ifndef bempp_elementary_linear_operator_hpp
#define bempp_elementary_linear_operator_hpp

#include "../common/common.hpp"

#include "linear_operator.hpp"

#include "../common/shared_ptr.hpp"

#include <stdexcept>
#include <vector>

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;
template <typename CoordinateType> class RawGridGeometry;
template <typename ValueType> class Basis;
class OpenClHandler;

} // namespace Fiber

namespace Bempp
{

/** \ingroup assembly
 *  \brief Linear operator whose weak form can be constructed with a single assembler.
 *
 *  The distinction between elementary and "non-elementary" linear operators is
 *  purely technical. An operator is considered elementary if its weak form can
 *  be constructed using a single instance of a subclass of
 *  Fiber::LocalAssemblerForOperators. Currently this is possible for the
 *  identity operator and for integral operators as defined in the
 *  documentation of ElementaryIntegralOperator.
 */
template <typename BasisFunctionType_, typename ResultType_>
class ElementaryLinearOperator : public LinearOperator<BasisFunctionType_, ResultType_>
{
    typedef LinearOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc LinearOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc LinearOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc LinearOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc LinearOperator::LocalAssemblerFactory */
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    /** \brief Type of the appropriate instantiation of Fiber::LocalAssemblerForOperators. */
    typedef Fiber::LocalAssemblerForOperators<ResultType> LocalAssembler;

    /** \copydoc LinearOperator::LinearOperator(const Space<BasisFunctionType>&, const Space<BasisFunctionType>&) */
    ElementaryLinearOperator(const Space<BasisFunctionType>& domain,
                             const Space<BasisFunctionType>& range,
                             const Space<BasisFunctionType>& dualToRange,
                             const std::string& label = "");

    ~ElementaryLinearOperator();

    virtual std::vector<const ElementaryLinearOperator<BasisFunctionType_, ResultType_>*>
    constituentOperators() const;
    virtual std::vector<ResultType_> constituentOperatorWeights() const;

    /** \brief Construct a local assembler suitable for this operator using a specified factory.
     *
     *  \param[in] assemblerFactory  Assembler factory to be used to construct the assembler.
     *
     *  (TODO: finish description of the other parameters.)
     */
    std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const;

    /** \brief Construct a local assembler suitable for this operator using a specified factory.
     *
     *  \param[in] assemblerFactory  Assembler factory to be used to construct the assembler.
     *  \param[in] options           Assembly options.
     *
     *  This is an overloaded function, provided for convenience. It
     *  automatically constructs most of the arguments required by the other
     *  overload.
     */
    std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const AssemblyOptions& options) const;

    /** \brief Assemble the operator's weak form using a specified local assembler.
     *
     *  This function is intended for internal use of the library. End users
     *  should not need to call it directly. They should use
     *  LinearOperator::assembleDetachedWeakForm() instead.
     */
    shared_ptr<DiscreteLinearOperator<ResultType_> >
    assembleWeakFormInternal(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry = NO_SYMMETRY);

private:
    virtual std::auto_ptr<LocalAssembler> makeAssemblerImpl(
            const LocalAssemblerFactory& assemblerFactory,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const = 0;

    virtual shared_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry) = 0;
};

} // namespace Bempp

#endif
