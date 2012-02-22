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

#ifndef fiber_standard_integration_manager_2d_hpp
#define fiber_standard_integration_manager_2d_hpp

#include "integration_manager.hpp"
#include "array_2d.hpp"
#include "element_pair_topology.hpp"
#include "opencl_options.hpp"
#include "separable_numerical_double_integrator.hpp"

#include "quadrature/galerkinduffy.hpp"
#include "quadrature/quadrature.hpp"

#include <boost/ptr_container/ptr_map.hpp>
#include <cstring>

namespace Fiber
{

template <typename ValueType, typename GeometryImp>
class StandardIntegrationManager2D :
        public IntegrationManager<ValueType, GeometryImp>
{    
private:
    struct DoubleIntegratorIndex
    {
        ElementPairTopology topology;
        int testOrder;
        int trialOrder;

        bool operator<(const DoubleIntegratorIndex& other) const {
            return topology < other.topology ?
                        true : testOrder < other.testOrder ?
                            true : trialOrder < other.trialOrder;
        }
    };

    typedef boost::ptr_map<DoubleIntegratorIndex,
    DoubleIntegrator<ValueType, GeometryImp> > IntegratorMap;

public:
    StandardIntegrationManager2D(
            const Expression<ValueType>& testExpression,
            const Kernel<ValueType>& kernel,
            const Expression<ValueType>& trialExpression,
            const OpenClOptions& openClOptions) :
        m_testExpression(testExpression),
        m_kernel(kernel),
        m_trialExpression(trialExpression),
        m_openClOptions(openClOptions)
    {}

    virtual const DoubleIntegrator<ValueType, GeometryImp>& testKernelTrialIntegrator(
            const GeometryImp& testGeometry,
            const GeometryImp& trialGeometry,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis) {
    if (testGeometry.dimension() != 2 || trialGeometry.dimension() != 2)
            throw std::invalid_argument("StandardIntegrationManager2D::"
                                        "testKernelTrialIntegrator: "
                                        "test and trial elements must be"
                                        "two-dimensional");

        DoubleIntegratorIndex integratorIndex;                
        integratorIndex.topology =
                determineElementPairTopology(testGeometry, trialGeometry);

        integratorIndex.testOrder = testBasis.order();
        integratorIndex.trialOrder = trialBasis.order();

        if (integratorIndex.topology.type = ElementPairTopology::Disjoint)
        {
            integratorIndex.testOrder +=
                    regularOrderIncrement(testGeometry, testBasis);
            integratorIndex.trialOrder +=
                    regularOrderIncrement(trialGeometry, trialBasis);
        }
        else // singular integral
        {
            integratorIndex.testOrder +=
                    singularOrderIncrement(testGeometry, testBasis);
            integratorIndex.trialOrder +=
                    singularOrderIncrement(trialGeometry, trialBasis);
        }

        return getIntegrator(integratorIndex);
    }

private:
    int regularOrderIncrement(const GeometryImp& geometry,
                               const Basis<ValueType>& basis) const
    {
        // Note: this function will make use of options supplied to the integrator
        // in its constructor

        // TODO:
        // 1. Check the size of elements and the distance between them
        //    and estimate the variability of the kernel
        // 2. Take into account the fact that elements might be isoparametric.

        // Make quadrature exact for a constant kernel and affine elements
        return 1;
    }

    int singularOrderIncrement(const GeometryImp& geometry,
                                const Basis<ValueType>& basis) const
    {
        // Note: this function will make use of options supplied to the integrator
        // in its constructor

        // TODO:
        // 1. Check the size of elements and estimate the variability of the
        //    (Sauter-Schwab-transformed) kernel
        // 2. Take into account the fact that elements might be isoparametric.

        // Make quadrature exact for a constant kernel and affine elements
        return 1;
    }

    const DoubleIntegrator<ValueType, GeometryImp>& getIntegrator(const DoubleIntegratorIndex& index)
    {
        const ElementPairTopology& topology = index.topology;

        typename IntegratorMap::iterator it = m_doubleIntegrators.find(index);
        if (it != m_doubleIntegrators.end())
            return *it->second;

        // Integrator doesn't exist yet and must be created.
        DoubleIntegrator<ValueType, GeometryImp>* integrator = 0;
        if (topology.type == ElementPairTopology::Disjoint)
        {
            // Create a tensor rule
            arma::Mat<ValueType> testPoints, trialPoints;
            std::vector<ValueType> testWeights, trialWeights;

            fillPointsAndWeightsRegular(topology.testVertexCount,
                                        index.testOrder,
                                        testPoints, testWeights);
            fillPointsAndWeightsRegular(topology.trialVertexCount,
                                        index.trialOrder,
                                        trialPoints, trialWeights);
            typedef SeparableNumericalDoubleIntegrator<ValueType, GeometryImp> Integrator;
            integrator = new Integrator(
                        testPoints, trialPoints, testWeights, trialWeights,
                        m_testExpression, m_kernel, m_trialExpression, m_openClOptions);
        }
        else
        {
            // We should create an appropriate non-tensor rule based on
            // Sauter-Schwab transformations and make a
            // NonseparableNumericalDoubleIntegrator, like this:

//            arma::Mat<ValueType> testPoints, trialPoints;
//            std::vector<ValueType> weights;

//            fillPointsAndWeightsSingular(index,
//                                         testPoints, trialPoints, weights);
//            integrator = new SingularDoubleIntegrator(
//                        testPoints, trialPoints, weights,
//                        m_testExpression, m_kernel, m_trialExpression,
//                        m_openClOptions);

            // However, in this first test version, we (incorrectly) create
            // a normal integrator.

            arma::Mat<ValueType> testPoints, trialPoints;
            std::vector<ValueType> testWeights, trialWeights;

            fillPointsAndWeightsRegular(topology.testVertexCount,
                                        index.testOrder,
                                        testPoints, testWeights);
            fillPointsAndWeightsRegular(topology.trialVertexCount,
                                        index.trialOrder,
                                        trialPoints, trialWeights);
            typedef SeparableNumericalDoubleIntegrator<ValueType, GeometryImp> Integrator;
            integrator = new Integrator(
                        testPoints, trialPoints, testWeights, trialWeights,
                        m_testExpression, m_kernel, m_trialExpression, m_openClOptions);

        }
        DoubleIntegratorIndex tmpIndex = index; // first argument of insert must be an l-value
        m_doubleIntegrators.insert(tmpIndex, integrator);
        return *integrator;
    }

    void fillPointsAndWeightsRegular(int vertexCount, int order,
                                     arma::Mat<ValueType>& points,
                                     std::vector<ValueType>& weights) const {
        if (vertexCount == 3)
            // Note that for triangles we want to have a different convention
            // for coordinates than Hyena does! TODO: implement this.
            reallyFillPointsAndWeightsRegular<TRIANGLE>(
                        order, points, weights);
        else
            return reallyFillPointsAndWeightsRegular<QUADRANGLE>(
                        order, points, weights);
    }

    template<ELEMENT_SHAPE SHAPE>
    void reallyFillPointsAndWeightsRegular(int order,
                                     arma::Mat<ValueType>& points,
                                     std::vector<ValueType>& weights) const {
        const int elementDim = 2;
        const QuadratureRule<SHAPE, GAUSS>& rule(order);
        const int pointCount = rule.getNumPoints();
        points.set_size(elementDim, pointCount);
        for (int i = 0; i < pointCount; ++i)
        {
            const Point2 point = rule.getPoint(i);
            for (int dim = 0; dim < elementDim; ++dim)
                points(dim, i) = point[dim];
            if (SHAPE == TRIANGLE)
                // Map points from
                // Hyena's reference triangle   (0,0)--(1,0)--(1,1)
                // to Dune's reference triangle (0,0)--(1,0)--(0,1).
                points(0, i) -= points(1, i);
        }
        weights.resize(pointCount);
        for (int i = 0; i < pointCount; ++i)
            weights[i] = rule.getWeight(i);
    }

// WORK IN PROGRESS
//    void fillPointsAndWeightsSingular(const DoubleIntegratorIndex& index,
//                                      Array2D<CoordinateType>& testPoints,
//                                      Array2D<CoordinateType>& trialPoints,
//                                      std::vector<ValueType>& weights) {
//        const ElementPairTopology& topology = index.topology;
//        if (topology.testVertexCount == 3 && topology.trialVertexCount == 3)
//        {
//            // Note that for triangles we want to have a different convention
//            // for coordinates than Hyena does! TODO: implement this.
//            if (topology.type == ElementPairTopology::SharedVertex)
//                reallyFillPointsAndWeightsSingular<TRIANGLE, VRTX_ADJACENT>(
//                            index, testPoints, trialPoints, weights);
//            else if (topology.type == ElementPairTopology::SharedEdge)
//                reallyFillPointsAndWeightsSingular<TRIANGLE, EDGE_ADJACENT>(
//                            index, testPoints, trialPoints, weights);
//            else if (topology.type == ElementPairTopology::Coincident)
//                reallyFillPointsAndWeightsSingular<TRIANGLE, COINCIDENT>(
//                            index, testPoints, trialPoints, weights);
//            else
//                throw std::invalid_argument("StandardIntegrationManager2D::"
//                                            "fillPointsAndWeightsSingular(): "
//                                            "Invalid element configuration. "
//                                            "This shouldn't happen!");
//        }
//        else if (topology.testVertexCount == 4 && topology.trialVertexCount == 4)
//        {
//            if (topology.type == ElementPairTopology::SharedVertex)
//                reallyFillPointsAndWeightsSingular<QUADRANGLE, VRTX_ADJACENT>(
//                            index, testPoints, trialPoints, weights);
//            else if (topology.type == ElementPairTopology::SharedEdge)
//                reallyFillPointsAndWeightsSingular<QUADRANGLE, EDGE_ADJACENT>(
//                            index, testPoints, trialPoints, weights);
//            else if (topology.type == ElementPairTopology::Coincident)
//                reallyFillPointsAndWeightsSingular<QUADRANGLE, COINCIDENT>(
//                            index, testPoints, trialPoints, weights);
//            else
//                throw std::invalid_argument("StandardIntegrationManager2D::"
//                                            "fillPointsAndWeightsSingular(): "
//                                            "Invalid element configuration. "
//                                            "This shouldn't happen!");
//        }
//        else
//            throw std::invalid_argument("StandardIntegrationManager2D::"
//                                        "fillPointsAndWeightsSingular(): "
//                                        "Singular quadrature rules for mixed "
//                                        "meshes are not implemented yet.");
//    }

//    template<ELEMENT_SHAPE SHAPE, SING_INT SINGULARITY>
//    void reallyFillPointsAndWeightsSingular(
//            const DoubleIntegratorIndex& index,
//            Array2D<CoordinateType>& testPoints,
//            Array2D<CoordinateType>& trialPoints,
//            std::vector<ValueType>& weights) {
//        // Is there a more efficient way?
//        const int order = std::max(index.testOrder, index.trialOrder);
//        const ElementaryQuadratureRule<SHAPE, GAUSS> rule(order);
//        const DuffyExpression<SHAPE, SINGULARITY> transform(rule);
//        const int pointCount = rule.getNumPoints();
//        const int regionCount = transform.getNumRegions();
//        const int totalPointCount = regionCount * pointCount;

//        Point2 point;

//        // Quadrature points
//        points.set_size(rule.shape_dim, totalPointCount);
//        for (int region = 0; region < regionCount; ++region)
//            for (int testIndex = 0; testIndex < pointCount; ++testIndex)
//                for (int trialIndex = 0; trialIndex < pointCount; ++trialIndex)
//                {
//                    int col = testIndex + trialIndex * pointCount +
//                            region * pointCount * pointCount;

//                    // Test point
//                    point = transform.getPointX(testIndex, trialIndex, region);
//                    for (int dim = 0; dim < rule.shape_dim; ++dim)
//                        testPoints(dim, col) = point.coords[dim];

//                    // Trial point
//                    point = transform.getPointY(testIndex, trialIndex, region);
//                    for (int dim = 0; dim < rule.shape_dim; ++dim)
//                        trialPoints(dim, col) = point.coords[dim];
//                }

//        // Weights
//        weights.resize(totalPointCount);
//        for (int region = 0; region < regionCount; ++region)
//            for (int testIndex = 0; testIndex < pointCount; ++testIndex)
//                for (int trialIndex = 0; trialIndex < pointCount; ++trialIndex)
//                    weights[i] = duffy.getWeight(testIndex, trialIndex, region) *
//                            rule.getWeight(testIndex) *
//                            rule.getWeight(trialIndex);

//        // TODO: "rotate" points if vertices are not in the appropriate positions.
//        if (SING_INT == VRTX_ADJACENT)
//        {
//            remapPointsSharedVertex<SHAPE>(index.topology.testSharedVertex0,
//                                           testPoints);
//            remapPointsSharedVertex<SHAPE>(index.topology.trialSharedVertex0,
//                                           trialPoints);
//        }
//        else if (SING_INT == EDGE_ADJACENT)
//        {
//            remapPointsSharedEdgeTestElement<SHAPE>(
//                        index.topology.testSharedVertex0,
//                        index.topology.testSharedVertex1,
//                        testPoints);
//            remapPointsSharedEdgeTrialElement<SHAPE>(
//                        index.topology.trialSharedVertex0,
//                        index.topology.trialSharedVertex1,
//                        trialPoints);
//        }
//    }

//    template <ELEMENT_SHAPE SHAPE> remapPointsSharedVertex(
//            int sharedVertex, Array2D<CoordinateType>& points) const
//    {} // default implementation

//    template <> remapPointsSharedVertex<TRIANGLE>(
//            int sharedVertex, Array2D<CoordinateType>& points) const
//    {
//        if (sharedVertex == 1)
//        {
//            const int pointCount = points.extents(1);
//            Array2D<CoordinateType> oldPoints(points);
//            for (int i = 0; i < pointCount; ++i)
//            {
//                points(0, i) = -oldPoints(0, i) - oldPoints(1, i) + 1.;
//                points(1, i) = oldPoints(1, i);
//            }
//        }
//        else if (sharedVertex == 2)
//        {
//            const int pointCount = points.extents(1);
//            Array2D<CoordinateType> oldPoints(points);
//            for (int i = 0; i < pointCount; ++i)
//            {
//                // unchecked
//                points(0, i) = oldPoints(0, i);
//                points(1, i) = -oldPoints(0, i) - oldPoints(1, i) + 1.;
//            }
//        }
//    }

//    template <ELEMENT_SHAPE SHAPE> remapPointsSharedEdgeTestElement(
//            int sharedVertex0, int sharedVertex1, Array2D<CoordinateType>& points) const
//    {} // default implementation , to be removed once everything's implemented

//    template <ELEMENT_SHAPE SHAPE> remapPointsSharedEdgeTrialElement(
//            int sharedVertex0, int sharedVertex1, Array2D<CoordinateType>& points) const
//    {} // default implementation , to be removed once everything's implemented

private:
    const Expression<ValueType>& m_testExpression;
    const Kernel<ValueType>& m_kernel;
    const Expression<ValueType>& m_trialExpression;
    OpenClOptions m_openClOptions;

    IntegratorMap m_doubleIntegrators;
};

} // namespace Fiber

#endif
