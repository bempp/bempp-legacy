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

#include "linear_operator.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/types.hpp"
#include "../fiber/types.hpp"

#include <vector>
#include <armadillo>

namespace Fiber
{

template <typename ValueType, typename GeometryImp> class IntegrationManager;
template <typename ValueType, typename GeometryImp> class IntegrationManagerFactory;

} // namespace Fiber

namespace Bempp
{

template <int codim> class EntityPointer;
class GeometryAdapter;

template <typename ValueType> class DiscreteScalarValuedLinearOperator;
template <typename ValueType> class DiscreteScalarValuedLinearOperator;
template <typename ValueType> class Space;
template <typename ValueType> class WeakFormAcaAssemblyHelper;

template <typename ValueType>
class ElementaryLinearOperator : public LinearOperator<ValueType>
{
    friend class WeakFormAcaAssemblyHelper<ValueType>;

public:
    typedef Fiber::IntegrationManagerFactory<ValueType, GeometryAdapter>
    IntegrationManagerFactory;
    typedef Fiber::IntegrationManager<ValueType, GeometryAdapter>
    IntegrationManager;

    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const IntegrationManagerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperator(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            const IntegrationManagerFactory& factory,
            const AssemblyOptions& options) const;

private:
    virtual std::auto_ptr<IntegrationManager >
    makeIntegrationManager(
            const IntegrationManagerFactory& factory) const = 0;

    /** \name Local assembly (virtual methods to be implemented
        in derived classes) @{ */
    /** \brief Assemble local weak forms.

    In this overload, a "column" of local weak forms is assembled. More
    specifically, on exit \p result is a vector of local weak forms corresponding
    to the following pairs of elements:

    - if \p callVariant is \p TEST_TRIAL, all pairs (\p elementA, \p elementB)
    for \p elementA in \p elementsA;

    - if \p callVariant is \p TRIAL_TEST, all pairs (\p elementB, \p elementA)
    for \p elementA in \p elementsA.

    Unless \p localDofIndexB is set to \p ALL_DOFS, only entries corresponding
    to the (\p localDofIndexB)th local DOF on \p elementB are calculated. */
    virtual void evaluateLocalWeakForms(
            CallVariant callVariant,
            const std::vector<const EntityPointer<0>*>& elementsA,
            const EntityPointer<0>& elementB,
            LocalDofIndex localDofIndexB,
            const Space<ValueType>& spaceA,
            const Space<ValueType>& spaceB,
            IntegrationManager& integrationMgr,
            std::vector<arma::Mat<ValueType> >& result) const = 0;

    /** \brief Assemble local weak forms.

    This overload constructs and assigns to the output parameter \p result the
    2D array of local weak forms corresponding to all pairs (testElement,
    trialElement) with testElement in testElements and trialElements in
    trialElements.

    This function should be used primarily for small blocks of elements lying
    close to each other. */
    virtual void evaluateLocalWeakForms(
            const std::vector<const EntityPointer<0>*>& testElements,
            const std::vector<const EntityPointer<0>*>& trialElements,
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            IntegrationManager& integrationMgr,
            Fiber::Array2D<arma::Mat<ValueType> >& result) const = 0;

    /** @}
        \name Operator assembly
        @{ */
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInDenseMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            IntegrationManager& integrationMgr,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInAcaMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            IntegrationManager& integrationMgr,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >

    /** @}
        \name Weak form assembly
        @{ */
    assembleWeakFormInDenseMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            IntegrationManager& integrationMgr,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInAcaMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            IntegrationManager& integrationMgr,
            const AssemblyOptions& options) const;
    /** @} */
};

} // namespace Bempp

#endif
