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

#ifndef bempp_fmm_far_field_helper_hpp
#define bempp_fmm_far_field_helper_hpp

#include <vector>
#include <complex>

#include "../common/armadillo_fwd.hpp"
#include "../common/common.hpp"
#include "../common/types.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
/** \endcond */

} // namespace Fiber


namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename ResultType> class Octree;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType> class LocalDofListsCache;
template <typename BasisFunctionType, typename ResultType> class FmmLocalAssembler;
template <typename ResultType> class FmmTransform;
class AssemblyOptions;
/** \endcond */

// fills the octree with data, depends on BasisFunctionType, whereas
// octree cannot, since it is needed by dicrete boundary operator

template <typename BasisFunctionType, typename ResultType>
class FmmFarFieldHelper
{
public:
    typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    FmmFarFieldHelper(
            const shared_ptr<Octree<ResultType> > octree,
            const Space<BasisFunctionType>& testSpace,
            const Space<BasisFunctionType>& trialSpace,
            const AssemblyOptions& options,
            const std::vector<unsigned int> &test_p2o,
            const std::vector<unsigned int> &trial_p2o,
            bool indexWithGlobalDofs,
            const FmmTransform<ResultType> &fmmTransform);

    void operator()(const tbb::blocked_range<unsigned int>& range) const;

private:
    arma::Mat<ResultType> evaluateFarFieldIntegrals(
        FmmLocalAssembler<BasisFunctionType, ResultType> &fmmLocalAssembler,
        const FmmTransform<ResultType> &fmmTransform,
        const arma::Col<CoordinateType> &nodeCentre, 
        const arma::Col<CoordinateType> &nodeSize, 
        unsigned int dofStartTrial, unsigned int dofCountTrial, bool isTest) const;

    shared_ptr<LocalDofListsCache<BasisFunctionType> > m_testDofListsCache;
    shared_ptr<LocalDofListsCache<BasisFunctionType> > m_trialDofListsCache;
    const AssemblyOptions& m_options;
    const Space<BasisFunctionType>& m_testSpace;
    const Space<BasisFunctionType>& m_trialSpace;

    shared_ptr<Octree<ResultType> > m_octree;
    const FmmTransform<ResultType> &m_fmmTransform;
};

} // namespace Bempp

#endif
