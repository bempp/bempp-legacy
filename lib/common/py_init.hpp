// Copyright (C) 2011-2015 by the BEM++ Authors
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


#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL bempp_ARRAY_API
#include "numpy/arrayobject.h"

#include <iostream>

#ifndef PY_INIT_HPP
#define PY_INIT_HPP

namespace Bempp {

class PyInit {

public:
  PyInit():
   m_pyFinalize(false) {
    if (!Py_IsInitialized()) {
      Py_Initialize();
      PyEval_InitThreads();
      m_pyFinalize = true;
      numpy_init();
    }
   }
  
  ~PyInit() {
      if (Py_IsInitialized() && m_pyFinalize) Py_Finalize();
  }

private:

  void* numpy_init() {
      import_array();
  }

  bool m_pyFinalize;
  PyInit(const PyInit &other);
  const PyInit &operator=(const PyInit &other);

  static PyInit m_singleton;
};
}

#endif
