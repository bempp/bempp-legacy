// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_standard_local_assembler_for_grid_functions_on_surfaces_hpp
#define fiber_standard_local_assembler_for_grid_functions_on_surfaces_hpp

#include "local_assembler_for_grid_functions.hpp"
#include "test_function_integrator.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <map>
#include <vector>

namespace Fiber
{

template <typename ValueType, typename IndexType> class OpenClHandler;

template <typename ValueType, typename GeometryFactory>
class StandardLocalAssemblerForGridFunctionsOnSurfaces :
        public LocalAssemblerForGridFunctions<ValueType>
{    
public:
    StandardLocalAssemblerForGridFunctionsOnSurfaces(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const Expression<ValueType>& testExpression,
            const Function<ValueType>& function,
            const OpenClHandler<ValueType,int>& openClHandler);
    virtual ~StandardLocalAssemblerForGridFunctionsOnSurfaces();

public:
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Col<ValueType> >& result);

private:
    const TestFunctionIntegrator<ValueType>& selectIntegrator(
            int elementIndex);

    int orderIncrement(int elementIndex) const;

    const TestFunctionIntegrator<ValueType>& getIntegrator(
            const SingleQuadratureDescriptor& index);

private:
    typedef tbb::concurrent_unordered_map<SingleQuadratureDescriptor,
    TestFunctionIntegrator<ValueType>*> IntegratorMap;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;
    const std::vector<const Basis<ValueType>*>& m_testBases;
    const Expression<ValueType>& m_testExpression;
    const Function<ValueType>& m_function;
    const OpenClHandler<ValueType,int>& m_openClHandler;

    IntegratorMap m_testFunctionIntegrators;
};

} // namespace Fiber

#include "standard_local_assembler_for_grid_functions_on_surfaces_imp.hpp"

#endif
