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

#ifndef bempp_auto_timer_hpp
#define bempp_auto_timer_hpp

#include <string>
#include <iostream>

#include <tbb/tick_count.h>

namespace Bempp {

/** \ingroup common
 *  \brief Timer that on destruction outputs the time elapsed since
 * construction. */
class AutoTimer {
public:
  /** \brief Constructor.

    \param[in] text Message to be printed on destruction. */
  explicit AutoTimer(const char *text = 0)
      : m_text(text), m_start(tbb::tick_count::now()) {}

  /** \overload */
  explicit AutoTimer(const std::string &text = std::string())
      : m_text(text), m_start(tbb::tick_count::now()) {}

  /** \brief Destructor. Print the previously specified message. */
  ~AutoTimer() {
    tbb::tick_count end = tbb::tick_count::now();
    std::cout << m_text << (end - m_start).seconds() << " s" << std::endl;
  }

private:
  std::string m_text;
  tbb::tick_count m_start;
};

} // namespace Bempp

#endif
