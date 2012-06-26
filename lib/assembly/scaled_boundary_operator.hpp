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

#ifndef bempp_scaled_boundary_operator_hpp
#define bempp_scaled_boundary_operator_hpp

#include "boundary_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename ResultType_>
class ScaledBoundaryOperator :
        public BoundaryOperator<BasisFunctionType_, ResultType_>
{
public:
    typedef BoundaryOperator<BasisFunctionType_, ResultType_> Base;
    /** \copydoc BoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc BoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc BoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc BoundaryOperator::LocalAssemblerFactory */
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;

    ScaledBoundaryOperator(ResultType weight, const Base& linOp);
    ScaledBoundaryOperator(const ScaledBoundaryOperator& other);

    virtual std::auto_ptr<BoundaryOperator<BasisFunctionType_, ResultType_> >
    clone() const;

    virtual bool supportsRepresentation(
            AssemblyOptions::Representation repr) const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                         const AssemblyOptions& options,
                         Symmetry symmetry);

private:
    ResultType m_weight;
    std::auto_ptr<Base> m_operator;
};

} // namespace Bempp

#endif
