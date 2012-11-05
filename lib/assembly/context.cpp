// Copyright (C) 2011 by the BEM++ Authors
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

#include "context.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/make_shared.hpp>
#include <iostream>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
        const shared_ptr<const QuadratureStrategy>& quadStrategy,
        const AssemblyOptions& assemblyOptions) :
    m_quadStrategy(quadStrategy),
    m_assemblyOptions(assemblyOptions)
{
    if (quadStrategy.get() == 0)
        throw std::invalid_argument("Context::Context(): "
                                    "quadStrategy must not be null");
}

template <typename BasisFunctionType, typename ResultType>
void Context<BasisFunctionType, ResultType>::printCachedOperatorStatistics() const
{
    std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType> > >
            cachedOperators = m_cache.aliveOperators();
    std::cout << cachedOperators.size() << " cached operators:\n";
    for (size_t i = 0; i < cachedOperators.size(); ++i)
        std::cout << i << " " << typeid(cachedOperators[i]).name() << std::endl;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp
