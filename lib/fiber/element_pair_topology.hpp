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

#ifndef fiber_element_pair_topology_hpp
#define fiber_element_pair_topology_hpp

#include <cassert>
#include <iostream>
#include <boost/tuple/tuple_comparison.hpp>

namespace Fiber
{

struct ElementPairTopology
{
    ElementPairTopology() :
        type(Disjoint),
        testVertexCount(0), trialVertexCount(0),
        testSharedVertex0(-1), testSharedVertex1(-1),
        trialSharedVertex0(-1), trialSharedVertex1(-1)
    {}

    enum Type { Disjoint, SharedVertex, SharedEdge, Coincident };
    Type type;
    unsigned char testVertexCount;
    unsigned char trialVertexCount;
    signed char testSharedVertex0;
    signed char testSharedVertex1;
    signed char trialSharedVertex0;
    signed char trialSharedVertex1;

    bool operator<(const ElementPairTopology& other) const {
        using boost::tuples::make_tuple;
        return make_tuple(type, testVertexCount, trialVertexCount,
                          testSharedVertex0, testSharedVertex1,
                          trialSharedVertex0, trialSharedVertex1) <
                make_tuple(other.type, other.testVertexCount, other.trialVertexCount,
                           other.testSharedVertex0, other.testSharedVertex1,
                           other.trialSharedVertex0, other.trialSharedVertex1);
    }

    bool operator==(const ElementPairTopology& other) const {
        return type == other.type &&
                testVertexCount == other.testVertexCount &&
                trialVertexCount == other.trialVertexCount &&
                testSharedVertex0 == other.testSharedVertex0 &&
                testSharedVertex1 == other.testSharedVertex1 &&
                trialSharedVertex0 == other.trialSharedVertex0 &&
                trialSharedVertex1 == other.trialSharedVertex1;
    }

    bool operator!=(const ElementPairTopology& other) const {
        return !operator==(other);
    }

    friend std::ostream&
    operator<< (std::ostream& dest, const ElementPairTopology& obj)
    {
        dest << obj.type << " "
             << (int)obj.testVertexCount << " "
             << (int)obj.trialVertexCount << " "
             << (int)obj.testSharedVertex0 << " "
             << (int)obj.testSharedVertex1 << " "
             << (int)obj.trialSharedVertex0 << " "
             << (int)obj.trialSharedVertex1;
        return dest;
    }
};

template <typename GeometryImp>
ElementPairTopology determineElementPairTopology(
        const GeometryImp& testGeometry, const GeometryImp& trialGeometry)
{
    ElementPairTopology topology;

    const int elementDim = testGeometry.dimension();
    assert(trialGeometry.dimension() == elementDim);

    // Retrieve indices of the vertices of test and trial elements
    const int MIN_VERTEX_COUNT = 2;
    const int MAX_VERTEX_COUNT = 4;
    topology.testVertexCount = testGeometry.vertexCount();
    topology.trialVertexCount = trialGeometry.vertexCount();
    if (topology.testVertexCount < MIN_VERTEX_COUNT ||
            MAX_VERTEX_COUNT < topology.testVertexCount ||
            topology.trialVertexCount < MIN_VERTEX_COUNT ||
            MAX_VERTEX_COUNT < topology.trialVertexCount)
        throw std::invalid_argument(
                "determineElementPairTopology: "
                "only linear, triangular and quadrilateral elements "
                "are supported");

//    GeometryImp::IndexType testVertexIndices[MAX_VERTEX_COUNT];
//    GeometryImp::IndexType trialVertexIndices[MAX_VERTEX_COUNT];
//    testGeometry.getVertexIndices(testVertexIndices);
//    trialGeometry.getVertexIndices(trialVertexIndices);

    // How many vertices coincide?
    int testSharedVertices[MAX_VERTEX_COUNT];
    int trialSharedVertices[MAX_VERTEX_COUNT];
    int hits = 0;

    for (int trialV = 0; trialV < topology.trialVertexCount; ++trialV)
        for (int testV = 0; testV < topology.testVertexCount; ++testV)
            if (testGeometry.vertexIndex(testV) ==
                    trialGeometry.vertexIndex(trialV))
            {
                testSharedVertices[hits] = testV;
                trialSharedVertices[hits] = trialV;
                ++hits;
                break;
            }

    if (hits == 0)
    {
        topology.type = ElementPairTopology::Disjoint;
    }
    else if (hits == 1)
    {
        topology.type = ElementPairTopology::SharedVertex;
        topology.testSharedVertex0 = testSharedVertices[0];
        topology.trialSharedVertex0 = trialSharedVertices[0];
    }
    else if (hits == 2 && elementDim == 2)
    {
        topology.type = ElementPairTopology::SharedEdge;
        topology.testSharedVertex0 = testSharedVertices[0];
        topology.testSharedVertex1 = testSharedVertices[1];
        topology.trialSharedVertex0 = trialSharedVertices[0];
        topology.trialSharedVertex1 = trialSharedVertices[1];
    }
    // coincident
    else if (hits == topology.testVertexCount &&
             hits == topology.trialVertexCount)
    {
        // Note: we don't handle the case where elements are different, but
        // their vertices coincide
        for (int i = 0; i < hits; ++i)
            assert(testSharedVertices[i] == trialSharedVertices[i]);
        topology.type = ElementPairTopology::Coincident;
    }
    else
        throw std::runtime_error(
                "Standard3DIntegrationManager::"
                "selectTestKernelTrialQuadratureRules() :"
                "Invalid element configuration");

    return topology;
}

} // namespace Fiber

#endif
