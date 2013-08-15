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

#ifndef bempp_weak_form_aca_assembly_helper_hpp
#define bempp_weak_form_aca_assembly_helper_hpp

#include "../common/common.hpp"

#include "ahmed_aux_fwd.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"
#include "../fiber/scalar_traits.hpp"

#include <tbb/atomic.h>
#include <tbb/concurrent_unordered_map.h>
#include <vector>

/** \cond FORWARD_DECL */
class cluster;
/** \endcond */

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
class AssemblyOptions;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class LocalDofListsCache;
template <typename BasisFunctionType> class Space;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief Class whose methods are called by Ahmed during assembly in the ACA mode.
 */
template <typename BasisFunctionType, typename ResultType>
class WeakFormAcaAssemblyHelper
{
public:
    typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;
    typedef Fiber::LocalAssemblerForOperators<ResultType> LocalAssembler;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef typename Fiber::ScalarTraits<ResultType>::RealType MagnitudeType;
    typedef typename AhmedTypeTraits<ResultType>::Type AhmedResultType;

    WeakFormAcaAssemblyHelper(const Space<BasisFunctionType>& testSpace,
                              const Space<BasisFunctionType>& trialSpace,
                              const std::vector<unsigned int>& p2oTestDofs,
                              const std::vector<unsigned int>& p2oTrialDofs,
                              const std::vector<LocalAssembler*>& assemblers,
                              const std::vector<const DiscreteLinOp*>& sparseTermsToAdd,
                              const std::vector<ResultType>& denseTermsMultipliers,
                              const std::vector<ResultType>& sparseTermsMultipliers,
                              const AssemblyOptions& options);

    /** \brief Evaluate entries of a general block.
     *
     *  Store the entries of the block defined
     *  by \p b1, \p n1, \p b2, \p n2 (in permuted ordering) in data. */
    void cmpbl(unsigned b1, unsigned n1, unsigned b2, unsigned n2,
               AhmedResultType* data,
               const cluster* c1 = 0, const cluster* c2 = 0,
               bool countAccessedEntries = true) const;

    /** \brief Evaluate entries of a symmetric block.
     *
     * Store the upper part of the (symmetric) block defined
     * by \p b1, \p n1, \p b1, \p n1 (in permuted ordering) columnwise in \p data. */
    void cmpblsym(unsigned b1, unsigned n1, AhmedResultType* data,
                  const cluster* c1 = 0,
                  bool countAccessedEntries = true) const;

    /** \brief Expected magnitude of the entries in this block. */
    MagnitudeType scale(unsigned b1, unsigned n1, unsigned b2, unsigned n2,
                        const cluster* c1 = 0, const cluster* c2 = 0) const;

    /** \brief Expected magnitude of the largest entry in this block relative to that
     *  of the largest entry in the whole matrix. */
    MagnitudeType relativeScale(unsigned b1, unsigned n1, unsigned b2, unsigned n2,
                                const cluster* c1 = 0, const cluster* c2 = 0) const;

    /** \brief Return the number of entries in the matrix that have been
     *  accessed so far. */
    size_t accessedEntryCount() const;

    /** \brief Reset the number of entries in the matrix that have been
     *  accessed so far. */
    void resetAccessedEntryCount();

private:
    MagnitudeType estimateMinimumDistance(
            const cluster* c1, const cluster* c2) const;

private:
    /** \cond PRIVATE */
    const Space<BasisFunctionType>& m_testSpace;
    const Space<BasisFunctionType>& m_trialSpace;
    const std::vector<unsigned int>& m_p2oTestDofs;
    const std::vector<unsigned int>& m_p2oTrialDofs;
    const std::vector<LocalAssembler*>& m_assemblers;
    const std::vector<const DiscreteLinOp*>& m_sparseTermsToAdd;
    const std::vector<ResultType>& m_denseTermsMultipliers;
    const std::vector<ResultType>& m_sparseTermsMultipliers;
    const AssemblyOptions& m_options;
    bool m_indexWithGlobalDofs;

    shared_ptr<LocalDofListsCache<BasisFunctionType> >
    m_testDofListsCache, m_trialDofListsCache;

    mutable tbb::atomic<size_t> m_accessedEntryCount;

    typedef tbb::concurrent_unordered_map<
    std::pair<const cluster*, const cluster*>, CoordinateType> DistanceMap;
    mutable DistanceMap m_distancesCache;

    /** \endcond */
};

} // namespace Bempp

#endif
