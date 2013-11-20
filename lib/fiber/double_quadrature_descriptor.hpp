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

#ifndef fiber_double_quadrature_descriptor_hpp
#define fiber_double_quadrature_descriptor_hpp

#include "element_pair_topology.hpp"

#include "../common/common.hpp"

#include <boost/tuple/tuple_comparison.hpp>
#include <ostream>

namespace Fiber
{

/** \brief Parameters of a quadrature rule used in the evaluation of
 *  integrals over pairs of elements. */
struct DoubleQuadratureDescriptor
{
    /** \brief Element pair configuration. */
    ElementPairTopology topology;
    /** \brief Degree of accuracy of the quadrature rule used on the test
     *  element. */
    int testOrder;
    /** \brief Degree of accuracy of the quadrature rule used on the trial
     *  element. */
    int trialOrder;

    bool operator<(const DoubleQuadratureDescriptor& other) const {
        using boost::tuples::make_tuple;
        return make_tuple(topology, testOrder, trialOrder) <
                make_tuple(other.topology, other.testOrder, other.trialOrder);
    }

    bool operator==(const DoubleQuadratureDescriptor& other) const {
        return topology == other.topology &&
                testOrder == other.testOrder &&
                trialOrder == other.trialOrder;
    }

    bool operator!=(const DoubleQuadratureDescriptor& other) const {
        return !operator==(other);
    }

    friend std::ostream&
    operator<< (std::ostream& dest, const DoubleQuadratureDescriptor& obj)
    {
        dest << obj.topology << " " << obj.testOrder << " " << obj.trialOrder;
        return dest;
    }
};

inline size_t tbb_hasher(const DoubleQuadratureDescriptor& d)
{
    const ElementPairTopology& t = d.topology;
    return (t.testVertexCount - 3) + 2 *
            ((t.trialVertexCount - 3) + 2 *
             (t.testSharedVertex0 + 4 *
              (t.trialSharedVertex0 + 4 *
               (t.testSharedVertex1 + 4 *
                (t.trialSharedVertex1 + 4 *
                 (d.testOrder + 256 *
                  d.trialOrder))))));
}

} // namespace Fiber

#endif
