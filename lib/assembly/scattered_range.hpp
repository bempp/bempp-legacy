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

#ifndef bempp_scattered_range_hpp
#define bempp_scattered_range_hpp

#include "../common/common.hpp"

#include <tbb/tbb_stddef.h> // split class
#include <cmath>
#include <cassert>

namespace Bempp {

class ScatteredRange {
public:
  class const_iterator {
  public:
    const_iterator(size_t value, size_t end, size_t step)
        : m_value(value), m_end(end), m_step(step) {
      assert(step > 0);
    }

    const_iterator operator++() {
      m_value = std::min(m_value + m_step, m_end);
      return *this;
    }

    size_t value() const { return m_value; }

    size_t step() const { return m_step; }

    // Implicit conversion to size_t
    operator size_t() const { return m_value; }

  private:
    size_t m_value, m_end, m_step;
  };

  ScatteredRange(size_t begin, size_t end, size_t step = 1)
      : m_begin(begin), m_end(end), m_step(step) {}

  ScatteredRange(ScatteredRange &range, tbb::split) {
    m_begin = range.m_begin + range.m_step;
    m_end = range.m_end;
    m_step = range.m_step * 2;
    range.m_step *= 2;
  }

  bool empty() const { return m_begin >= m_end; }

  bool is_divisible() const { return m_begin + m_step < m_end; }

  const_iterator begin() const {
    return const_iterator(m_begin, m_end, m_step);
  }

  const_iterator end() const { return const_iterator(m_end, m_end, m_step); }

  size_t step() const { return m_step; }

  size_t size() const {
    if (empty())
      return 0;
    else
      return (m_end - 1 - m_begin) / m_step + 1;
  }

private:
  size_t m_begin, m_end, m_step;
};

} // namespace Bempp

#endif
