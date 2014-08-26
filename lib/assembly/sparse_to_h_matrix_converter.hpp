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

#ifndef bempp_sparse_to_h_matrix_converter_hpp
#define bempp_sparse_to_h_matrix_converter_hpp

#include "../common/common.hpp"
#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "ahmed_aux_fwd.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/boost_shared_array_fwd.hpp"

#include <vector>

namespace Bempp {

// This structure is intended for internal BEM++ use only.
template <typename ValueType> struct SparseToHMatrixConverter {
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
  typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
  typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
  typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

  static void constructHMatrix(int *rowOffsets, int *colIndices, double *values,
                               std::vector<unsigned int> &domain_o2p,
                               std::vector<unsigned int> &range_p2o, double eps,
                               AhmedBemBlcluster *blockCluster,
                               boost::shared_array<AhmedMblock *> &mblocks,
                               int &maximumRank);
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
