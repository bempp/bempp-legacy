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

#ifndef bempp_lazy_hpp
#define bempp_lazy_hpp

// After Intel(R) Threading Building Blocks Design Patterns, modified

#include <tbb/atomic.h>
#include <tbb/mutex.h>
#include <memory>

namespace Bempp
{

/** \brief Thread-safe wrapper of a lazily inititialised object.

  \tparam T           Type of the stored object.
  \tparam Initialiser Type of a copy-constructible functor providing the
                      member function
                        std::auto_ptr<T> operator()()
                      that will be used to construct the object the first time
                      it is requested.
  \tparam Mutex       Type of the mutex to be used to lock the object during
                      initialisation.
*/
template<typename T, typename Initialiser, typename Mutex=tbb::mutex>
class Lazy
{
public:
    Lazy(const Initialiser& init) :
        m_init(init), m_value() {
    }

    ~Lazy() {
        delete m_value;
    }

    T& get() {
        if (!m_value) {
            typename Mutex::scoped_lock lock(m_mutex);
            if (!m_value) {
                std::auto_ptr<T> t = m_init();
                m_value = t.release();
            }
        }
        return *m_value;
    }

private:
    Initialiser m_init;
    tbb::atomic<T*> m_value;
    Mutex m_mutex;
};

} // namespace Bempp

#endif
