// Copyright (C) 2011-2014 by the BEM++ Authors
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

#include "weak_form_hmat_assembly_helper.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType>::
    WeakFormHMatAssemblyHelper(
        const Space<BasisFunctionType> &testSpace,
        const Space<BasisFunctionType> &trialSpace,
        const shared_ptr<hmat::BlockClusterTree<>> blockClusterTree,
        const std::vector<LocalAssembler *> &assemblers,
        const std::vector<const DiscreteLinOp *> &sparseTermsToAdd,
        const std::vector<ResultType> &denseTermsMultipliers,
        const std::vector<ResultType> &sparseTermsMultipliers,
        const AssemblyOptions &options)
    : m_testSpace(testSpace), m_trialSpace(trialSpace),
      m_blockClusterTree(blockClusterTree), m_assemblers(assemblers),
      m_sparseTermsToAdd(sparseTermsToAdd),
      m_denseTermsMultipliers(denseTermsMultipliers),
      m_sparseTermsMultipliers(sparseTermsMultipliers), m_options(options) {}
}
