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

#ifndef bempp_linear_operator_superposition_hpp
#define bempp_linear_operator_superposition_hpp

#include "linear_operator.hpp"

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class ElementaryLinearOperator;

/** \brief Superposition of linear operators.
 *
 *  \ingroup assembly
 */
template <typename BasisFunctionType, typename ResultType>
class LinearOperatorSuperposition :
        public LinearOperator<BasisFunctionType, ResultType>
{
public:
    typedef LinearOperator<BasisFunctionType, ResultType> Base;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef typename Fiber::LocalAssemblerForOperators<ResultType>
    LocalAssembler;

    LinearOperatorSuperposition(const Base& term1, const Base& term2);

    LinearOperatorSuperposition(const Base& term, const ResultType& scalar);

    virtual int trialComponentCount() const;
    virtual int testComponentCount() const;

    virtual std::vector<const ElementaryLinearOperator<BasisFunctionType, ResultType>*>
    constituentOperators() const;
    virtual std::vector<ResultType> constituentOperatorWeights() const;

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

private:
    /** @name Constituent elementary operators list management
     *  @{ */

    /** \brief Append operators to the list of constituent elementary operators.
     *
     *  \param[in] operators
     *    Vector of pointers to the elementary linear operators to be appended
     *    to the list of constituent operators. These objects must continue to
     *    exist at least until the weak form of this operator is assembled.
     *
     *  \param[in] weights
     *    Vector of the corresponding weights.
     *
     *  \see constituentOperators(), constituentOperatorWeights().
     */
    void addConstituentOperators(
            const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>&
            operators,
            const std::vector<ResultType>& weights);

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormImpl(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInDenseMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInAcaMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInArbitraryMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

private:
    std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
    m_constituentOperators;
    std::vector<ResultType> m_constituentOperatorWeights;
};

} //namespace Bempp

#endif
