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

#include "serial_blas_region.hpp"

#include "bempp/common/config_blas_and_lapack.hpp"

#if defined(WITH_MKL)

// #include <mkl.h>
extern "C" {
void MKL_Set_Num_Threads(int nth);
int MKL_Get_Max_Threads(void);
#define mkl_set_num_threads MKL_Set_Num_Threads
#define mkl_get_max_threads MKL_Get_Max_Threads
} // extern "C"

// Currently, dynamic changing of the number of threads used by
// GotoBLAS/OpenBLAS
// does not work; the user needs to set the environmental variable
// GOTO_NUM_THREADS to 1 before running programs linked to BEM++
//#elif defined(WITH_GOTOBLAS) || defined(WITH_OPENBLAS)
// extern "C" {
//#include <common.h>
//}
#endif

#include <iostream>

namespace Fiber {

SerialBlasRegion::SerialBlasRegion() {
#if defined(WITH_MKL)
  m_originalThreadCount = mkl_get_max_threads();
  mkl_set_num_threads(1);
//#elif defined(WITH_GOTOBLAS) || defined(WITH_OPENBLAS)
//    m_originalThreadCount = get_num_procs();
//    goto_set_num_threads(1);
#else
  m_originalThreadCount = -1;
#endif
}

SerialBlasRegion::~SerialBlasRegion() {
#if defined(WITH_MKL)
  mkl_set_num_threads(m_originalThreadCount);
//#elif defined(WITH_GOTOBLAS) || defined(WITH_OPENBLAS)
//    goto_set_num_threads(m_originalThreadCount);
#endif
}

} // namespace Fiber
