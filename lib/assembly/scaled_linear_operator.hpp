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

#ifndef bempp_scaled_linear_operator_hpp
#define bempp_scaled_linear_operator_hpp

#include "linear_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename ResultType_>
class ScaledLinearOperator :
        public LinearOperator<BasisFunctionType_, ResultType_>
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

    ScaledLinearOperator(ResultType weight, const Base& linOp);
    ScaledLinearOperator(const ScaledLinearOperator& other);

    virtual std::auto_ptr<LinearOperator<BasisFunctionType_, ResultType_> >
    clone() const;

    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType_> >
    assembleDetachedWeakFormImpl(const LocalAssemblerFactory& factory,
                                 const AssemblyOptions& options,
                                 Symmetry symmetry) const;

private:
    ResultType m_weight;
    std::auto_ptr<Base> m_operator;
};

} // namespace Bempp

#endif
