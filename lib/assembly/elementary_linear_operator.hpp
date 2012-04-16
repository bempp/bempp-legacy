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

#include <vector>
#include <stdexcept>

namespace Fiber
{

template <typename ValueType> class LocalAssemblerForOperators;
template <typename ValueType> class RawGridGeometry;
template <typename ValueType> class Basis;
template <typename CoordinateType, typename IndexType> class OpenClHandler;

} // namespace Fiber

namespace Bempp
{

template <typename ValueType>
class ElementaryLinearOperator : public LinearOperator<ValueType>
{
public:
    typedef typename LinearOperator<ValueType>::LocalAssemblerFactory
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForOperators<ValueType> LocalAssembler;

    ElementaryLinearOperator(const Space<ValueType>& testSpace, const Space<ValueType>& trialSpace)
        : LinearOperator<ValueType>(testSpace,trialSpace){
        std::vector<ElementaryLinearOperator<ValueType> const*> v;
        std::vector<ValueType> m;
        v.push_back(this);
        m.push_back(1.0);
        addLocalOperatorsMultipliers(v,m);

    }

    /** \brief Using a specified factory, construct a local assembler suitable
       for this operator. */
    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const GeometryFactory& geometryFactory,
            const Fiber::RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Fiber::Basis<ValueType>*>& testBases,
            const std::vector<const Fiber::Basis<ValueType>*>& trialBases,
            const Fiber::OpenClHandler<ValueType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** \brief Assemble the operator's weak form using a specified local assembler.

      This function is not intended to be called directly by the user. */
    virtual std::auto_ptr<DiscreteLinearOperator<ValueType> >
    assembleWeakFormInternal(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const = 0;

//    /** \brief Multiply the operator in-place by a scalar.

//      This method affects the results of subsequent calls to assembleWeakForm()
//      and assembleOperator(). */
//    void scale(ValueType multiplier) {
//        m_multiplier = multiplier;
//    }

//    /** \brief Return the current value of the scalar by which this operator is multiplied. */
//    ValueType multiplier() const {
//        return m_multiplier;
//    }


};

}

#endif // bempp_elementary_linear_operator_hpp
