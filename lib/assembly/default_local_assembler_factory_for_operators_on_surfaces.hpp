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

#ifndef bempp_default_local_assembler_factory_for_operators_on_surfaces_hpp
#define bempp_default_local_assembler_factory_for_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "../fiber/default_local_assembler_factory_for_operators_on_surfaces.hpp"
#include "../grid/geometry_factory.hpp"

namespace Bempp
{

using Fiber::AccuracyOptions;

template <typename BasisFunctionType, typename ResultType>
class DefaultLocalAssemblerFactoryForOperatorsOnSurfaces :
        public Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<
        BasisFunctionType, ResultType, GeometryFactory>
{
private:
    typedef Fiber::DefaultLocalAssemblerFactoryForOperatorsOnSurfaces<
    BasisFunctionType, ResultType, GeometryFactory> Base;
public:
    /** \brief Construct a local assembler factory with default accuracy settings. */
    DefaultLocalAssemblerFactoryForOperatorsOnSurfaces();

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit DefaultLocalAssemblerFactoryForOperatorsOnSurfaces(
            const AccuracyOptions& accuracyOptions);
};

} // namespace Bempp

#endif
