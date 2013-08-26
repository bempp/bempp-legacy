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

#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "weak_form_aca_assembly_helper.hpp"

#include "assembly_options.hpp"
#include "discrete_boundary_operator.hpp"
#include "local_dof_lists_cache.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../grid/grid_view.hpp"
#include "../space/space.hpp"

#include "ahmed_aux.hpp"
#include "ahmed_complex.hpp"

#include <map>
#include <set>
#include <utility>
#include <cstdio>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::WeakFormAcaAssemblyHelper(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<unsigned int>& p2oTestDofs,
        const std::vector<unsigned int>& p2oTrialDofs,
        const std::vector<LocalAssembler*>& assemblers,
        const std::vector<const DiscreteLinOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermsMultipliers,
        const std::vector<ResultType>& sparseTermsMultipliers,
        const AssemblyOptions& options) :
    m_testSpace(testSpace), m_trialSpace(trialSpace),
    m_p2oTestDofs(p2oTestDofs), m_p2oTrialDofs(p2oTrialDofs),
    m_assemblers(assemblers), m_sparseTermsToAdd(sparseTermsToAdd),
    m_denseTermsMultipliers(denseTermsMultipliers),
    m_sparseTermsMultipliers(sparseTermsMultipliers),
    m_options(options),
    m_indexWithGlobalDofs(m_options.acaOptions().mode != AcaOptions::HYBRID_ASSEMBLY),
    m_testDofListsCache(new LocalDofListsCache<BasisFunctionType>(
                            m_testSpace, m_p2oTestDofs, m_indexWithGlobalDofs)),
    m_trialDofListsCache(new LocalDofListsCache<BasisFunctionType>(
                            m_trialSpace, m_p2oTrialDofs, m_indexWithGlobalDofs))
    //,
//    m_trialDofListsCache(&testSpace == &trialSpace &&
//                         std::equal(p2oTestDofs.begin(), p2oTestDofs.end(),
//                                    p2oTrialDofs.begin()) ?
//                             m_testDofListsCache :
//                             boost::make_shared<LocalDofListsCache<BasisFunctionType> >(
//                                 m_trialSpace, m_p2oTrialDofs, m_indexWithGlobalDofs))
{
    if (!m_indexWithGlobalDofs && !m_sparseTermsToAdd.empty())
        throw std::invalid_argument(
                "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
                "combining sparse and dense terms in hybrid ACA mode "
                "is not supported at present");
    for (size_t i = 0; i < assemblers.size(); ++i)
        if (!assemblers[i])
            throw std::invalid_argument(
                    "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
                    "no elements of the 'assemblers' vector may be null");
    for (size_t i = 0; i < sparseTermsToAdd.size(); ++i)
        if (!sparseTermsToAdd[i])
            throw std::invalid_argument(
                    "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
                    "no elements of the 'sparseTermsToAdd' vector may be null");
    resetAccessedEntryCount();
}

template <typename BasisFunctionType, typename ResultType>
typename WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::estimateMinimumDistance(
        const cluster* c1, const cluster* c2) const
{
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;

    AhmedBemCluster* cluster1 =
        const_cast<AhmedBemCluster*>(// AHMED is not const-correct
            dynamic_cast<const AhmedBemCluster*>(c1));
    AhmedBemCluster* cluster2 =
        const_cast<AhmedBemCluster*>(// AHMED is not const-correct
            dynamic_cast<const AhmedBemCluster*>(c2));
    // Lower bound on the minimum distance between elements from the two clusters
    CoordinateType minDist = -1.; // negative, read: unknown

   if (cluster1 && cluster2)
       minDist = sqrt(cluster1->extDist2(cluster2));

    // else
        // std::cout << "Warning: clusters not available" << std::endl;




    // if (cluster1 && cluster2) {
    //     std::pair<const cluster*, const cluster*> key(cluster1, cluster2);
    //     typename DistanceMap::const_iterator it = m_distancesCache.find(key);
    //     if (it != m_distancesCache.end())
    //         minDist = it->second;
    //     else {
    //         shared_ptr<const LocalDofLists<BasisFunctionType> > testDofLists =
    //                 m_testDofListsCache->get(c1->getnbeg(), c1->size());
    //         shared_ptr<const LocalDofLists<BasisFunctionType> > trialDofLists =
    //                 m_trialDofListsCache->get(c2->getnbeg(), c2->size());
    //         minDist = m_assemblers[0]->estimateMinimumDistance(
    //                     testDofLists->elementIndices, trialDofLists->elementIndices);
    //         m_distancesCache.insert(std::make_pair(key, minDist));
    //     }
    // }

    return minDist;
}

template <typename BasisFunctionType, typename ResultType>
void WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::cmpbl(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        AhmedResultType* ahmedData,
        const cluster* c1, const cluster* c2, bool countAccessedEntries) const
{
    //    std::cout << "\nRequested block: (" << b1 << ", " << n1 << "; "
    //              << b2 << ", " << n2 << ")" << std::endl;
    if (countAccessedEntries)
        m_accessedEntryCount += n1 * n2;

    // if negative, it means: unknown
    const CoordinateType minDist = estimateMinimumDistance(c1, c2);

    // This is a non-op for real types. For complex types, it converts a pointer
    // to Ahmed's scomp (resp. dcomp) to a pointer to std::complex<float>
    // (resp. std::complex<double>), which should be perfectly safe since these
    // types have the same binary representation.
    ResultType* data = reinterpret_cast<ResultType*>(ahmedData);

    // Convert AHMED matrix indices into DOF indices
    shared_ptr<const LocalDofLists<BasisFunctionType> > testDofLists =
            m_testDofListsCache->get(b1, n1);
    shared_ptr<const LocalDofLists<BasisFunctionType> > trialDofLists =
            m_trialDofListsCache->get(b2, n2);

    // Requested original matrix indices
    typedef typename LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
    const std::vector<DofIndex>& testOriginalIndices =
            testDofLists->originalIndices;
    const std::vector<DofIndex>& trialOriginalIndices =
            trialDofLists->originalIndices;
    // Necessary elements
    const std::vector<int>& testElementIndices =
            testDofLists->elementIndices;
    const std::vector<int>& trialElementIndices =
            trialDofLists->elementIndices;
    // Necessary local dof indices in each element
    const std::vector<std::vector<LocalDofIndex> >& testLocalDofs =
            testDofLists->localDofIndices;
    const std::vector<std::vector<LocalDofIndex> >& trialLocalDofs =
            trialDofLists->localDofIndices;
    // Weights of local dofs in each element
    const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights =
            testDofLists->localDofWeights;
    const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights =
            trialDofLists->localDofWeights;
    for (size_t i = 0; i < testLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < testLocalDofWeights[i].size(); ++j)
            assert(std::abs(testLocalDofWeights[i][j]) > 0.);
    for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
            assert(std::abs(trialLocalDofWeights[i][j]) > 0.);

    // Corresponding row and column indices in the matrix to be calculated
    // and stored in ahmedData
    const std::vector<std::vector<int> >& blockRows =
            testDofLists->arrayIndices;
    const std::vector<std::vector<int> >& blockCols =
            trialDofLists->arrayIndices;

    arma::Mat<ResultType> result(data, n1, n2, false /*copy_aux_mem*/,
                                 true /*strict*/);
    result.fill(0.);

    // First, evaluate the contributions of the dense terms
    if (n2 == 1)
    {
        // Only one column of the block needed. This means that we need only
        // one local DOF from just one or a few trialElements. Evaluate the
        // local weak form for one local trial DOF at a time.

        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t nTrialElem = 0;
             nTrialElem < trialElementIndices.size();
             ++nTrialElem)
        {
            const int activeTrialElementIndex = trialElementIndices[nTrialElem];

            // The body of this loop will very probably only run once (single
            // local DOF per trial element)
            for (size_t nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof)
            {
                LocalDofIndex activeTrialLocalDof =
                        trialLocalDofs[nTrialElem][nTrialDof];
                BasisFunctionType activeTrialLocalDofWeight =
                        trialLocalDofWeights[nTrialElem][nTrialDof];
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
                {
                    m_assemblers[nTerm]->evaluateLocalWeakForms(
                                Fiber::TEST_TRIAL, testElementIndices,
                                activeTrialElementIndex, activeTrialLocalDof,
                                localResult, minDist);
                    for (size_t nTestElem = 0;
                         nTestElem < testElementIndices.size();
                         ++nTestElem)
                        for (size_t nTestDof = 0;
                             nTestDof < testLocalDofs[nTestElem].size();
                             ++nTestDof)
                            result(blockRows[nTestElem][nTestDof], 0) +=
                                    m_denseTermsMultipliers[nTerm] *
                                    conj(testLocalDofWeights[nTestElem][nTestDof]) *
                                    activeTrialLocalDofWeight *
                                    localResult[nTestElem](testLocalDofs[nTestElem][nTestDof]);
                }
            }
        }
    }
    else if (n1 == 1) // very few testElements
    {
        // Only one row of the block needed. This means that we need only
        // one local DOF from just one or a few testElements. Evaluate the
        // local weak form for one local test DOF at a time.

        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t nTestElem = 0;
             nTestElem < testElementIndices.size();
             ++nTestElem)
        {
            const int activeTestElementIndex = testElementIndices[nTestElem];
            // The body of this loop will very probably only run once (single
            // local DOF per test element)
            for (size_t nTestDof = 0;
                 nTestDof < testLocalDofs[nTestElem].size();
                 ++nTestDof)
            {
                LocalDofIndex activeTestLocalDof =
                        testLocalDofs[nTestElem][nTestDof];
                BasisFunctionType activeTestLocalDofWeight =
                        testLocalDofWeights[nTestElem][nTestDof];
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
                {
                    m_assemblers[nTerm]->evaluateLocalWeakForms(
                                Fiber::TRIAL_TEST, trialElementIndices,
                                activeTestElementIndex, activeTestLocalDof,
                                localResult, minDist);
                    for (size_t nTrialElem = 0;
                         nTrialElem < trialElementIndices.size();
                         ++nTrialElem)
                        for (size_t nTrialDof = 0;
                             nTrialDof < trialLocalDofs[nTrialElem].size();
                             ++nTrialDof)
                            result(0, blockCols[nTrialElem][nTrialDof]) +=
                                    m_denseTermsMultipliers[nTerm] *
                                    conj(activeTestLocalDofWeight) *
                                    trialLocalDofWeights[nTrialElem][nTrialDof] *
                                    localResult[nTrialElem](trialLocalDofs[nTrialElem][nTrialDof]);
                }
            }
        }
    }
    else if (n1 <= 32 && n2 <= 32) // a "fat" block
    {
        // The whole block or its submatrix needed. This means that we are
        // likely to need all or almost all local DOFs from most elements.
        // Evaluate the full local weak form for each pair of test and trial
        // elements and then select the entries that we need.

        Fiber::_2dArray<arma::Mat<ResultType> > localResult;
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
        {
            m_assemblers[nTerm]->evaluateLocalWeakForms(
                    testElementIndices, trialElementIndices, localResult,
                    minDist);
            for (size_t nTrialElem = 0;
                 nTrialElem < trialElementIndices.size();
                 ++nTrialElem)
                for (size_t nTrialDof = 0;
                     nTrialDof < trialLocalDofs[nTrialElem].size();
                     ++nTrialDof)
                    for (size_t nTestElem = 0;
                         nTestElem < testElementIndices.size();
                         ++nTestElem)
                        for (size_t nTestDof = 0;
                             nTestDof < testLocalDofs[nTestElem].size();
                             ++nTestDof)
                            result(blockRows[nTestElem][nTestDof],
                                   blockCols[nTrialElem][nTrialDof]) +=
                                    m_denseTermsMultipliers[nTerm] *
                                    conj(testLocalDofWeights[nTestElem][nTestDof]) *
                                    trialLocalDofWeights[nTrialElem][nTrialDof] *
                                    localResult(nTestElem, nTrialElem)
                                    (testLocalDofs[nTestElem][nTestDof],
                                     trialLocalDofs[nTrialElem][nTrialDof]);
        }
    }
    else
    {
        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t nTestElem = 0;
             nTestElem < testElementIndices.size();
             ++nTestElem)
        {
            const int activeTestElementIndex = testElementIndices[nTestElem];
            // The body of this loop will very probably only run once (single
            // local DOF per test element)
            for (size_t nTestDof = 0;
                 nTestDof < testLocalDofs[nTestElem].size();
                 ++nTestDof)
            {
                LocalDofIndex activeTestLocalDof =
                        testLocalDofs[nTestElem][nTestDof];
                BasisFunctionType activeTestLocalDofWeight =
                        testLocalDofWeights[nTestElem][nTestDof];
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
                {
                    m_assemblers[nTerm]->evaluateLocalWeakForms(
                                Fiber::TRIAL_TEST, trialElementIndices,
                                activeTestElementIndex, activeTestLocalDof,
                                localResult, minDist);
                    for (size_t nTrialElem = 0;
                         nTrialElem < trialElementIndices.size();
                         ++nTrialElem)
                        for (size_t nTrialDof = 0;
                             nTrialDof < trialLocalDofs[nTrialElem].size();
                             ++nTrialDof)
                            result(blockRows[nTestElem][nTestDof],
                                   blockCols[nTrialElem][nTrialDof]) +=
                                    m_denseTermsMultipliers[nTerm] *
                                    conj(activeTestLocalDofWeight) *
                                    trialLocalDofWeights[nTrialElem][nTrialDof] *
                                    localResult[nTrialElem](trialLocalDofs[nTrialElem][nTrialDof]);
                }
            }
        }
    }

    // Probably can be removed
    if (m_indexWithGlobalDofs)
        // Now, add the contributions of the sparse terms
        for (size_t nTerm = 0; nTerm < m_sparseTermsToAdd.size(); ++nTerm)
            m_sparseTermsToAdd[nTerm]->addBlock(
                        // since m_indexWithGlobalDofs is set, these refer
                        // to global DOFs
                        testOriginalIndices, trialOriginalIndices,
                        m_sparseTermsMultipliers[nTerm], result);
    // else m_sparseTermsToAdd is empty (as we have verified in the constructor)
}

template <typename BasisFunctionType, typename ResultType>
void WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::cmpblsym(
        unsigned b1, unsigned n1, AhmedResultType* ahmedData,
        const cluster* c1, bool countAccessedEntries) const
{
    if (countAccessedEntries)
        m_accessedEntryCount += n1 * n1;

    // Calculate results as usual (without taking symmetry into account)
    arma::Mat<ResultType> block(n1, n1);
    cmpbl(b1, n1, b1, n1, ahmedCast(block.memptr()), c1, c1);

    // and now store the upper part of the matrix in the memory block
    // provided by Ahmed
    for (size_t col = 0; col < n1; ++col)
        for (size_t row = 0; row <= col; ++row)
            *ahmedData++ = ahmedCast(block(row, col));
}

template <typename BasisFunctionType, typename ResultType>
typename WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::scale(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        const cluster* c1, const cluster* c2) const
{
    typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;

    AhmedBemCluster* cluster1 =
        const_cast<AhmedBemCluster*>(// AHMED is not const-correct
            dynamic_cast<const AhmedBemCluster*>(c1));
    AhmedBemCluster* cluster2 =
        const_cast<AhmedBemCluster*>(// AHMED is not const-correct
            dynamic_cast<const AhmedBemCluster*>(c2));
    if (cluster1 && cluster2) {
        // both getdiam2() and dist2() are effectively const, but not declared so
        // in AHMED
        const CoordinateType dist = sqrt(cluster1->extDist2(cluster2));
        MagnitudeType result = 0.;
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
            result = std::max(result,
                              m_assemblers[nTerm]->estimateRelativeScale(dist));
        return result;
    }
    else
        return m_options.acaOptions().scaling;
}

template <typename BasisFunctionType, typename ResultType>
typename WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::relativeScale(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        const cluster* c1, const cluster* c2) const
{
    const CoordinateType minDist = estimateMinimumDistance(c1, c2);
    if (minDist < 0)
        return 1.;
    MagnitudeType result = 0.;
    for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
        result = std::max(result,
                          m_assemblers[nTerm]->estimateRelativeScale(minDist));
    return result;
}

template <typename BasisFunctionType, typename ResultType>
size_t
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::
accessedEntryCount() const
{
    return m_accessedEntryCount;
}

template <typename BasisFunctionType, typename ResultType>
void
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::
resetAccessedEntryCount()
{
    m_accessedEntryCount = 0;
}

// Explicit instantiations

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(WeakFormAcaAssemblyHelper);

}

#endif // WITH_AHMED
