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

template <typename ArgumentType, typename ResultType>
class ElementaryLinearOperator;

// only scalar multipliers allowed, tensor ones would
// require knowledge of vector components distribution
// in the discrete operator
template <typename ArgumentType, typename ResultType>
class LinearOperatorSuperposition :
        public LinearOperator<ArgumentType, ResultType>
{
public:
    typedef LinearOperator<ArgumentType, ResultType> Base;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef typename Fiber::LocalAssemblerForOperators<ResultType>
    LocalAssembler;

    LinearOperatorSuperposition(const Base& term1, const Base& term2);

    LinearOperatorSuperposition(const Base& term, const ResultType& scalar);

    virtual int trialComponentCount() const;
    virtual int testComponentCount() const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakForm(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

private:
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInDenseMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInAcaMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInArbitraryMode(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;
};

} //namespace Bempp

#endif
