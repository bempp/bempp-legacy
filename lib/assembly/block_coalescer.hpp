// Copyright (C) 2013 by the BEM++ Authors
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

#ifndef bempp_block_coalescer_hpp
#define bempp_block_coalescer_hpp

#include "../common/common.hpp"
#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "aca_options.hpp"
#include "ahmed_aux_fwd.hpp"
#include "../common/boost_scoped_array_fwd.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

class Epetra_CrsMatrix;

namespace Bempp {

template <typename ValueType> class BlockCoalescer {
  typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
  typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
  typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
  typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

public:
  BlockCoalescer(blcluster *blclustersRoot, blcluster *decomposedBlclustersRoot,
                 const shared_ptr<const Epetra_CrsMatrix> &
                     permutedTestGlobalToFlatLocalMap,
                 const shared_ptr<const Epetra_CrsMatrix> &
                     permutedTrialGlobalToFlatLocalMap,
                 const boost::shared_array<AhmedMblock *> &blocks,
                 const boost::shared_array<AhmedMblock *> &decomposedBlocks,
                 const AcaOptions &acaOptions);

  void coalesceBlock(unsigned index);

private:
  /** \cond PRIVATE */
  void coalesceDenseBlock(unsigned index);
  void coalesceLowRankBlock(unsigned index);

private:
  boost::scoped_array<blcluster *> m_blclusters;
  boost::scoped_array<blcluster *> m_decomposedBlclusters;
  shared_ptr<const Epetra_CrsMatrix> m_permutedTestGlobalToFlatLocalMap;
  shared_ptr<const Epetra_CrsMatrix> m_permutedTrialGlobalToFlatLocalMap;
  boost::shared_array<AhmedMblock *> m_blocks;
  boost::shared_array<AhmedMblock *> m_decomposedBlocks;
  AcaOptions m_acaOptions;
  /** \endcond PRIVATE */
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
