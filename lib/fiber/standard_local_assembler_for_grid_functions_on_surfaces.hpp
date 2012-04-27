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

template <typename BasisFunctionType, typename FunctionValueType, typename GeometryFactory>
class StandardLocalAssemblerForGridFunctionsOnSurfaces :
        public LocalAssemblerForGridFunctions<
        typename Coercion<BasisFunctionType, FunctionValueType>::Type>
{    
public:
    typedef typename Coercion<BasisFunctionType, FunctionValueType>::Type ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    StandardLocalAssemblerForGridFunctionsOnSurfaces(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const Expression<CoordinateType>& testExpression,
            const Function<FunctionValueType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler);
    virtual ~StandardLocalAssemblerForGridFunctionsOnSurfaces();

public:
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Col<ResultType> >& result);

private:
    typedef TestFunctionIntegrator<BasisFunctionType, FunctionValueType> Integrator;

    const Integrator& selectIntegrator(int elementIndex);

    int orderIncrement(int elementIndex) const;

    const Integrator& getIntegrator(const SingleQuadratureDescriptor& index);

private:
    typedef tbb::concurrent_unordered_map<SingleQuadratureDescriptor,
    Integrator*> IntegratorMap;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<CoordinateType>& m_rawGeometry;
    const std::vector<const Basis<BasisFunctionType>*>& m_testBases;
    const Expression<CoordinateType>& m_testExpression;
    const Function<FunctionValueType>& m_function;
    const OpenClHandler<CoordinateType, int>& m_openClHandler;

    IntegratorMap m_testFunctionIntegrators;
};

} // namespace Fiber

#endif
