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

#ifndef fiber_quadrature_options_hpp
#define fiber_quadrature_options_hpp

#include "../common/common.hpp"

namespace Fiber
{

/** \brief Options controlling the order of numerical quadrature.
 *
 *  This class can be used to specify the required order of accuracy of a
 *  numerical quadrature rule used to approximate a certain class of integrals.
 *  The order can be set either as absolute or as relative with respect to some
 *  default value determined by the code doing the integration. */
class QuadratureOptions
{
public:
    /** \brief Construct a QuadratureOptions object corresponding
     *  to a default accuracy order. */
    QuadratureOptions() : m_relative(true), m_value(0) {
    }

    /** \brief Construct a QuadratureOptions object corresponding
     *  to a specific accuracy order. */
    QuadratureOptions(int order, bool relative) :
        m_relative(relative), m_value(order) {
    }

    /** \brief Set quadrature accuracy order to \p order. */
    void setAbsoluteQuadratureOrder(int order) {
        m_relative = false;
        m_value = order;
    }

    /** \brief Set quadrature accuracy order to a default value plus \p offset. */
    void setRelativeQuadratureOrder(int offset) {
        m_relative = true;
        m_value = offset;
    }

    /** \brief Get quadrature accuracy order assuming that its default value is
     *  \p defaultOrder. */
    int quadratureOrder(int defaultOrder) const {
        if (m_relative)
            return defaultOrder + m_value;
        else
            return m_value;
    }

private:
    bool m_relative;
    int m_value;
};

} // namespace Fiber

#endif
