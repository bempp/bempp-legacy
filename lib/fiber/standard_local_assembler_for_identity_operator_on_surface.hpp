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

#ifndef fiber_standard_local_assembler_for_identity_operator_on_surface_hpp
#define fiber_standard_local_assembler_for_identity_operator_on_surface_hpp

#include "local_assembler_for_identity_operator.hpp"
#include "opencl_options.hpp"
#include "numerical_quadrature.hpp"
#include "numerical_test_trial_integrator.hpp"

#include <armadillo>
#include <boost/static_assert.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
class StandardLocalAssemblerForIdentityOperatorOnSurface :
        public LocalAssemblerForIdentityOperator<ValueType>
{    
public:
    StandardLocalAssemblerForIdentityOperatorOnSurface(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<ValueType>& rawGeometry,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            const Expression<ValueType>& testExpression,
            const Expression<ValueType>& trialExpression,
            const OpenClOptions& openClOptions) :
        m_geometryFactory(geometryFactory),
        m_rawGeometry(rawGeometry),
        m_testBases(testBases),
        m_trialBases(trialBases),
        m_testExpression(testExpression),
        m_trialExpression(trialExpression),
        m_openClOptions(openClOptions)
    {
        if (rawGeometry.vertices().n_rows != 3)
            throw std::invalid_argument(
                    "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                    "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                    "vertex coordinates must be three-dimensional");
        const int elementCount = rawGeometry.elementCornerIndices().n_cols;
        if (rawGeometry.elementCornerIndices().n_rows < 3 ||
                4 < rawGeometry.elementCornerIndices().n_rows)
            throw std::invalid_argument(
                    "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                    "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                    "Elements must have either 3 or 4 corners");
        if (!rawGeometry.auxData().is_empty() &&
                rawGeometry.auxData().n_cols != elementCount)
            throw std::invalid_argument(
                    "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                    "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                    "number of columns of auxData must match that of "
                    "elementCornerIndices");
        if (testBases.size() != elementCount)
            throw std::invalid_argument(
                    "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                    "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                    "size of testBases must match the number of columns of "
                    "elementCornerIndices");
        if (trialBases.size() != elementCount)
            throw std::invalid_argument(
                    "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                    "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                    "size of trialBases must match the number of columns of "
                    "elementCornerIndices");
    }

public:
    virtual void evaluateLocalWeakForms(
            const std::vector<int>& elementIndices,
            std::vector<arma::Mat<ValueType> >& result);

private:
    const TestTrialIntegrator<ValueType>& selectIntegrator(int elementIndex) {
        SingleQuadratureDescriptor desc;

        // Get number of corners of the specified element
        desc.vertexCount = m_rawGeometry.elementCornerCount(elementIndex);

        // Determine integrand's order and required quadrature order
        const int expressionOrder =
                m_testBases[elementIndex]->order() +
                m_trialBases[elementIndex]->order();
        desc.order = ((expressionOrder + 1) + 1 /* round up */) / 2;

        return getIntegrator(desc);
    }

private:
    const TestTrialIntegrator<ValueType>& getIntegrator(
            const SingleQuadratureDescriptor& desc)
    {
        typename IntegratorMap::iterator it = m_testTrialIntegrators.find(desc);
        if (it != m_testTrialIntegrators.end())
        {
//            std::cout << "getIntegrator(: " << index << "): integrator found" << std::endl;
            return *it->second;
        }
//        std::cout << "getIntegrator(: " << index << "): integrator not found" << std::endl;

        // Integrator doesn't exist yet and must be created.
        arma::Mat<ValueType> points;
        std::vector<ValueType> weights;
        fillSingleQuadraturePointsAndWeights(desc.vertexCount, desc.order,
                                             points, weights);

        typedef NumericalTestTrialIntegrator<ValueType, GeometryFactory> Integrator;
        std::auto_ptr<TestTrialIntegrator<ValueType> > integrator(
                new Integrator(points, weights,
                               m_geometryFactory, m_rawGeometry,
                               m_testExpression, m_trialExpression,
                               m_openClOptions));

        return *m_testTrialIntegrators.insert(desc, integrator).first->second;
    }

private:
    typedef boost::ptr_map<SingleQuadratureDescriptor,
    TestTrialIntegrator<ValueType> > IntegratorMap;

private:
    const GeometryFactory& m_geometryFactory;
    const RawGridGeometry<ValueType>& m_rawGeometry;
    const std::vector<const Basis<ValueType>*>& m_testBases;
    const std::vector<const Basis<ValueType>*>& m_trialBases;
    const Expression<ValueType>& m_testExpression;
    const Expression<ValueType>& m_trialExpression;
    OpenClOptions m_openClOptions;

    IntegratorMap m_testTrialIntegrators;
};

} // namespace Fiber

#include "standard_local_assembler_for_identity_operator_on_surface_imp.hpp"

#endif
