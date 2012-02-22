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
class Geometry;

template <typename ValueType> class DiscreteScalarValuedLinearOperator;
template <typename ValueType> class DiscreteScalarValuedLinearOperator;
template <typename ValueType> class Space;
template <typename ValueType> class WeakFormAcaAssemblyHelper;

template <typename ValueType>
class ElementaryLinearOperator : public LinearOperator<ValueType>
{
    friend class WeakFormAcaAssemblyHelper<ValueType>;

public:
    virtual std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakForm(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const Fiber::IntegrationManagerFactory<ValueType, Geometry>& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperator(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            const Fiber::IntegrationManagerFactory<ValueType, Geometry>& factory,
            const AssemblyOptions& options) const;

private:
    virtual std::auto_ptr<Fiber::IntegrationManager<ValueType, Geometry> > makeIntegrationManager(
            const Fiber::IntegrationManagerFactory<ValueType, Geometry>& factory) const = 0;

    /** \name Local assembly (virtual methods to be implemented
        in derived classes) @{ */
    virtual void evaluateLocalWeakForms(
            CallVariant evalVariant,
            const std::vector<const EntityPointer<0>*>& elementsA,
            const EntityPointer<0>& elementB,
            LocalDofIndex localDofIndexB,
            const Space<ValueType>& spaceA,
            const Space<ValueType>& spaceB,
            Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
            std::vector<arma::Mat<ValueType> >& result) const = 0;

    /** @}
        \name Operator assembly
        @{ */
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInDenseMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
    assembleOperatorInAcaMode(
            const arma::Mat<ctype>& testPoints,
            const Space<ValueType>& trialSpace,
            Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
            const AssemblyOptions& options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >

    /** @}
        \name Weak form assembly
        @{ */
    assembleWeakFormInDenseMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
    assembleWeakFormInAcaMode(
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
            const AssemblyOptions& options) const;
    /** @} */
};

} // namespace Bempp

#endif
