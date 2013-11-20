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

#ifndef fiber_single_quadrature_descriptor_hpp
#define fiber_single_quadrature_descriptor_hpp

#include "../common/common.hpp"

#include <boost/tuple/tuple_comparison.hpp>
#include <ostream>

namespace Fiber
{

/** \brief Parameters of a quadrature rule used in the evaluation of
 *  integrals over single elements. */
struct SingleQuadratureDescriptor
{
    /** \brief Number of vertices of the element constituing the integration domain. */
    int vertexCount;
    /** \brief Degree of accuracy of the quadrature rule. */
    int order;

    bool operator<(const SingleQuadratureDescriptor& other) const {
        return std::make_pair(vertexCount, order) <
                std::make_pair(other.vertexCount, other.order);
    }

    bool operator==(const SingleQuadratureDescriptor& other) const {
        return vertexCount == other.vertexCount &&
                order == other.order;
    }

    bool operator!=(const SingleQuadratureDescriptor& other) const {
        return !operator==(other);
    }

    friend std::ostream&
    operator<< (std::ostream& dest, const SingleQuadratureDescriptor& obj)
    {
        dest << obj.vertexCount << " " << obj.order;
        return dest;
    }
};

inline size_t tbb_hasher(const SingleQuadratureDescriptor& d)
{
    return (d.vertexCount - 3) + 2 * d.order;
}

} // namespace Fiber

#endif
