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

#ifndef bempp_source_term_hpp
#define bempp_source_term_hpp

#include <memory>

namespace Fiber
{

template <typename ValueType> class Expression;
template <typename ValueType, typename GeometryFactory>
class LocalAssemblerFactory;
template <typename ValueType> class Function;
template <typename ValueType> class LocalAssemblerForSourceTerms;

} // namespace Fiber

namespace Bempp
{

class AssemblyOptions;
class GeometryFactory;
template <typename ValueType> class DiscreteScalarValuedSourceTerm;
template <typename ValueType> class Function;
template <typename ValueType> class Space;

template <typename ValueType>
class SourceTerm
{
public:
    typedef Fiber::LocalAssemblerFactory<ValueType, GeometryFactory>
    LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForSourceTerms<ValueType> LocalAssembler;

    ~SourceTerm() {}

    std::auto_ptr<DiscreteScalarValuedSourceTerm<ValueType> >
    assembleWeakForm(
            const Fiber::Function<ValueType>& function,
            const Space<ValueType>& testSpace,
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

private:
    std::auto_ptr<DiscreteScalarValuedSourceTerm<ValueType> >
    reallyAssembleWeakForm(
            const Space<ValueType>& testSpace,
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
};

} // namespace Bempp

#endif
