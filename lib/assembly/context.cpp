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
//#include "numerical_quadrature_strategy.hpp"

#include <boost/make_shared.hpp>

namespace Bempp
{

//template <typename BasisFunctionType, typename ResultType>
//shared_ptr<const Context<BasisFunctionType, ResultType> >
//Context<BasisFunctionType, ResultType>::m_defaultContext(
//        boost::make_shared<Context<BasisFunctionType, ResultType>(
//            AssemblyOptions(),
//            boost::make_shared<DefaultQuadratureStrategyForOperatorsOnSurfaces<
//            BasisFunctionType, ResultType> >()));

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
        const shared_ptr<QuadratureStrategy>& quadStrategy,
        const AssemblyOptions& assemblyOptions) :
    m_quadStrategy(quadStrategy),
    m_assemblyOptions(assemblyOptions)
{
    if (quadStrategy.get() == 0)
        throw std::invalid_argument("Context::Context(): "
                                    "quadStrategy cannot be null");
}

//template <typename BasisFunctionType, typename ResultType>
//void Context<BasisFunctionType, ResultType>::setDefaultContext(
//        const shared_ptr<const Context>& newContext)
//{
//    if (newContext.get() == 0)
//        throw std::invalid_argument("Context::setDefaultContext(): "
//                                    "new default context cannot be null");
//    m_defaultContext = newContext;
//}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp
