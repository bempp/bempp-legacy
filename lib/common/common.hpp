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

#ifndef bempp_common_hpp
#define bempp_common_hpp

// Needs to be included before Python support
#include <boost/property_tree/ptree.hpp>

#ifdef __INTEL_COMPILER
#pragma warning(disable : 279 186 858 262 1011 654 597)
#endif

#include <Python.h>
// Yes, use the new syntax from Dune 2.2
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1
#include "bempp/common/bempp_dune_config.hpp"

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL bempp_ARRAY_API
#include "numpy/arrayobject.h"

#include <stdexcept>

#endif
