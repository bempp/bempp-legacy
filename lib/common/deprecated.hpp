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

#ifndef bempp_deprecated_hpp
#define bempp_deprecated_hpp

#include "common.hpp"

/** \ingroup common
 *  \def BEMPP_DEPRECATED
 *  \brief Macro used to mark deprecated functions or classes.
 *
 *  Functions and classes marked as BEMPP_DEPRECATED are liable to be removed
 *  in future versions of the library. */
#if defined(__llvm__)
#define BEMPP_DEPRECATED
#else
#if defined(__GNUC__)
#define BEMPP_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define BEMPP_DEPRECATED __declspec(deprecated)
#else
#define BEMPP_DEPRECATED
#endif
#endif

// Macros for temporarily disabling deprecation warnings, by Jonathan Wakely
// (http://gcc.gnu.org/ml/gcc-help/2011-01/msg00135.html); modified.
// Should be called like
// BEMPP_GCC_DIAG_OFF(deprecated)
// BEMPP_GCC_DIAG_ON(deprecated)
#if defined(__llvm__)
#define BEMPP_GCC_DIAG_OFF(x)
#define BEMPP_GCC_DIAG_ON(x)
#else
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
#define BEMPP_GCC_DIAG_STR(s) #s
#define BEMPP_GCC_DIAG_JOINSTR(x, y) BEMPP_GCC_DIAG_STR(x##y)
#define BEMPP_GCC_DIAG_DO_PRAGMA(x) _Pragma(#x)
#define BEMPP_GCC_DIAG_PRAGMA(x) BEMPP_GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#define BEMPP_GCC_DIAG_OFF(x)                                                  \
  BEMPP_GCC_DIAG_PRAGMA(push)                                                  \
      BEMPP_GCC_DIAG_PRAGMA(ignored BEMPP_GCC_DIAG_JOINSTR(-W, x))
#define BEMPP_GCC_DIAG_ON(x) BEMPP_GCC_DIAG_PRAGMA(pop)
#else
#define BEMPP_GCC_DIAG_OFF(x)                                                  \
  BEMPP_GCC_DIAG_PRAGMA(ignored BEMPP_GCC_DIAG_JOINSTR(-W, x))
// This doesn't seem to work for GCC 4.4, so we just give up and don't restore
// the warning level.
// #  define BEMPP_GCC_DIAG_ON(x)  BEMPP_GCC_DIAG_PRAGMA(warning
// BEMPP_GCC_DIAG_JOINSTR(-W,x))
#define BEMPP_GCC_DIAG_ON(x)
#endif
#else
#define BEMPP_GCC_DIAG_OFF(x)
#define BEMPP_GCC_DIAG_ON(x)
#endif
#endif

#endif // bempp_deprecated_hpp
