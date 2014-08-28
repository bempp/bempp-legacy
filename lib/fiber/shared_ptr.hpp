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

#ifndef fiber_shared_ptr_hpp
#define fiber_shared_ptr_hpp

#include "../common/common.hpp"
#include <functional>

#ifdef __INTEL_COMPILER
#pragma warning(disable : 279 858)
#endif

#include <boost/shared_ptr.hpp>
#include <boost/pointer_cast.hpp>

#ifdef __INTEL_COMPILER
#pragma warning(default : 858)
#endif

// Extend std::hash to deal with boost shared pointers

namespace std {

template <typename T> struct std::hash<boost::shared_ptr<T>> {
  std::size_t operator()(const boost::shared_ptr<T> &x) const {
    return std::hash<T *>()(x.get());
  }
};
}

namespace Fiber {

using boost::shared_ptr;
using boost::dynamic_pointer_cast;

struct null_deleter {
  void operator()(const void *) const {}
};

/** \brief Create a shared pointer from a reference to an object allocated on
stack.

The object will not be deleted when the shared pointer goes out of scope. */
template <typename T> inline shared_ptr<T> make_shared_from_ref(T &t) {
  return shared_ptr<T>(&t, null_deleter());
}

template <typename T>
inline shared_ptr<const T> make_shared_from_const_ref(const T &t) {
  return shared_ptr<const T>(&t, null_deleter());
}

} // namespace Fiber

#endif // fiber_shared_ptr_hpp
