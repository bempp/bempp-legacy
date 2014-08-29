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

#ifndef bempp_hmat_parameters_hpp
#define bempp_hmat_parameters_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterEntry.hpp>

namespace Bempp {

using Teuchos::ParameterList;
using Teuchos::ParameterEntry;

/** \ingroup parameters
 *  \brief Parameters controlling the HMat library.
 */

struct HMatParameters {

  /** \brief Specifies whether to assemble using global or
   *  local dofs. */
  static ParameterEntry defaultHMatAssemblyMode();

  /** \brief Specifies the minimum block size below which
   *  blocks are assumed to be dense.
   *  \note Default Value: (unsigned int) 16 */
  static ParameterEntry defaultMinBlockSize();

  /** \brief Specifies the maximum size of an admissible block.
   * \note Default Value: (unsigned int) 2048 */
  static ParameterEntry defaultMaxBlockSize();

  /** \brief Specifies the block separation parameter eta.
   * \note Default Value: (double) 1.2 */
  static ParameterEntry defaultEta();

  /** \brief Create a \p ParameterList object containing the
   *  default options. */
  static shared_ptr<ParameterList> parameterList();
};
}
#endif
