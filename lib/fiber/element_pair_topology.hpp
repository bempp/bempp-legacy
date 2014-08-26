// Copyright (C) 2011-2012 by the Bem++ Authors
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

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include <cassert>
#include <iostream>
#include <boost/tuple/tuple_comparison.hpp>

namespace Fiber {

/** \brief Configuration of a pair of elements. */
struct ElementPairTopology {
  /** \brief Constructor. */
  ElementPairTopology()
      : type(Disjoint), testVertexCount(0), trialVertexCount(0),
        testSharedVertex0(-1), testSharedVertex1(-1), trialSharedVertex0(-1),
        trialSharedVertex1(-1) {}

  /** \brief Location of one element with respect to the other. */
  enum Type {
    Disjoint,
    SharedVertex,
    SharedEdge,
    Coincident
  };
  /** \brief Location of one element with respect to the other. */
  Type type;
  /** \brief Number of vertices of the test element. */
  unsigned char testVertexCount;
  /** \brief Number of vertices of the trial element. */
  unsigned char trialVertexCount;
  /** \brief Index of the first vertex of the test element that is
   *  shared with the trial element, or -1 if no vertices are
   *  shared. */
  signed char testSharedVertex0;
  /** \brief Index of the second vertex of the test element that is
   *  shared with the trial element, or -1 if at most one vertex is
   *  shared. */
  signed char testSharedVertex1;
  /** \brief Index of the first vertex of the trial element that is
   *  shared with the test element, or -1 if no vertices are
   *  shared. */
  signed char trialSharedVertex0;
  /** \brief Index of the second vertex of the trial element that is
   *  shared with the test element, or -1 if at most one vertex is
   *  shared. */
  signed char trialSharedVertex1;

  bool operator<(const ElementPairTopology &other) const {
    using boost::tuples::make_tuple;
    return make_tuple(type, testVertexCount, trialVertexCount,
                      testSharedVertex0, testSharedVertex1, trialSharedVertex0,
                      trialSharedVertex1) <
           make_tuple(other.type, other.testVertexCount, other.trialVertexCount,
                      other.testSharedVertex0, other.testSharedVertex1,
                      other.trialSharedVertex0, other.trialSharedVertex1);
  }

  bool operator==(const ElementPairTopology &other) const {
    return type == other.type && testVertexCount == other.testVertexCount &&
           trialVertexCount == other.trialVertexCount &&
           testSharedVertex0 == other.testSharedVertex0 &&
           testSharedVertex1 == other.testSharedVertex1 &&
           trialSharedVertex0 == other.trialSharedVertex0 &&
           trialSharedVertex1 == other.trialSharedVertex1;
  }

  bool operator!=(const ElementPairTopology &other) const {
    return !operator==(other);
  }

  friend std::ostream &operator<<(std::ostream &dest,
                                  const ElementPairTopology &obj) {
    dest << obj.type << " " << (int)obj.testVertexCount << " "
         << (int)obj.trialVertexCount << " " << (int)obj.testSharedVertex0
         << " " << (int)obj.testSharedVertex1 << " "
         << (int)obj.trialSharedVertex0 << " " << (int)obj.trialSharedVertex1;
    return dest;
  }
};

inline ElementPairTopology determineElementPairTopologyIn3D(
    const arma::Col<int> &testElementCornerIndices,
    const arma::Col<int> &trialElementCornerIndices) {
  ElementPairTopology topology;

// Determine number of element corners
#ifndef NDEBUG
  const int MIN_VERTEX_COUNT = 3;
#endif
  const int MAX_VERTEX_COUNT = 4;
  topology.testVertexCount = testElementCornerIndices.n_rows;
  assert(MIN_VERTEX_COUNT <= topology.testVertexCount &&
         topology.testVertexCount <= MAX_VERTEX_COUNT);
  topology.trialVertexCount = trialElementCornerIndices.n_rows;
  assert(MIN_VERTEX_COUNT <= topology.trialVertexCount &&
         topology.trialVertexCount <= MAX_VERTEX_COUNT);

  // How many vertices coincide?
  int testSharedVertices[MAX_VERTEX_COUNT];
  int trialSharedVertices[MAX_VERTEX_COUNT];
  int hits = 0;

  for (int trialV = 0; trialV < topology.trialVertexCount; ++trialV)
    for (int testV = 0; testV < topology.testVertexCount; ++testV)
      if (testElementCornerIndices(testV) ==
          trialElementCornerIndices(trialV)) {
        testSharedVertices[hits] = testV;
        trialSharedVertices[hits] = trialV;
        ++hits;
        break;
      }

  if (hits == 0) {
    topology.type = ElementPairTopology::Disjoint;
  } else if (hits == 1) {
    topology.type = ElementPairTopology::SharedVertex;
    topology.testSharedVertex0 = testSharedVertices[0];
    topology.trialSharedVertex0 = trialSharedVertices[0];
  } else if (hits == 2) // && elementDim == 2)
  {
    topology.type = ElementPairTopology::SharedEdge;
    topology.testSharedVertex0 = testSharedVertices[0];
    topology.testSharedVertex1 = testSharedVertices[1];
    topology.trialSharedVertex0 = trialSharedVertices[0];
    topology.trialSharedVertex1 = trialSharedVertices[1];
  }
  // coincident
  else if (hits == topology.testVertexCount &&
           hits == topology.trialVertexCount) {
    // Note: we don't handle the case where elements are different, but
    // their vertices coincide
    for (int i = 0; i < hits; ++i)
      assert(testSharedVertices[i] == trialSharedVertices[i]);
    topology.type = ElementPairTopology::Coincident;
  } else
    throw std::runtime_error("Standard3DIntegrationManager::"
                             "selectTestKernelTrialQuadratureRules() :"
                             "Invalid element configuration");

  return topology;
}

} // namespace Fiber

#endif
