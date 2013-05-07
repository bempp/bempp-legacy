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

#include "potential_operator_aca_assembly_helper.hpp"

#include "evaluation_options.hpp"
#include "discrete_boundary_operator.hpp"
#include "local_dof_lists_cache.hpp"
#include "component_lists_cache.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../fiber/types.hpp"
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
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::
PotentialOperatorAcaAssemblyHelper(
        const arma::Mat<CoordinateType>& points,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<unsigned int>& p2oPoints,
        const std::vector<unsigned int>& p2oTrialDofs,
        const std::vector<LocalAssembler*>& assemblers,
        const std::vector<ResultType>& termMultipliers,
        const EvaluationOptions& options) :
    m_points(points), m_trialSpace(trialSpace),
    m_p2oPoints(p2oPoints), m_p2oTrialDofs(p2oTrialDofs),
    m_assemblers(assemblers),
    m_termMultipliers(termMultipliers),
    m_options(options),
    m_indexWithGlobalDofs(m_options.acaOptions().globalAssemblyBeforeCompression),
    m_trialDofListsCache(new LocalDofListsCache<BasisFunctionType>(
                            m_trialSpace, m_p2oTrialDofs, m_indexWithGlobalDofs))
{
    if (assemblers.empty())
        throw std::invalid_argument(
                "PotentialOperatorAcaAssemblyHelper::"
                "PotentialOperatorAcaAssemblyHelper(): "
                "the 'assemblers' vector must not be null");
    if (assemblers.size() != termMultipliers.size())
        throw std::invalid_argument(
                "PotentialOperatorAcaAssemblyHelper::"
                "PotentialOperatorAcaAssemblyHelper(): "
                "the 'assemblers' and 'termMultipliers' vectors must have the "
                "same length");
    for (size_t i = 0; i < assemblers.size(); ++i)
        if (!assemblers[i])
            throw std::invalid_argument(
                    "PotentialOperatorAcaAssemblyHelper::"
                    "PotentialOperatorAcaAssemblyHelper(): "
                    "no elements of the 'assemblers' vector may be null");
    m_componentCount = assemblers[0]->resultDimension();
    for (size_t i = 1; i < assemblers.size(); ++i)
        if (assemblers[i]->resultDimension() != m_componentCount)
            throw std::invalid_argument(
                    "PotentialOperatorAcaAssemblyHelper::"
                    "PotentialOperatorAcaAssemblyHelper(): "
                    "all assemblers must produce results with the same number "
                    "of components");
    m_componentListsCache.reset(
                new ComponentListsCache(m_p2oPoints, m_componentCount));
    resetAccessedEntryCount();
}

template <typename BasisFunctionType, typename ResultType>
typename PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::estimateMinimumDistance(
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
    if (cluster1 && cluster2) {
        // both getdiam2() and dist2() are effectively const, but not declared so
        // in AHMED
        const double diam1 = sqrt(cluster1->getdiam2());
        const double diam2 = sqrt(cluster2->getdiam2());
        const double dist = sqrt(cluster1->extDist2(cluster2));
        minDist = std::max(0., dist - (diam1 + diam2) / 2.);
    }
    // else
        // std::cout << "Warning: clusters not available" << std::endl;
    return minDist;
}

template <typename BasisFunctionType, typename ResultType>
void PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::cmpbl(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        AhmedResultType* ahmedData,
        const cluster* c1, const cluster* c2, bool countAccessedEntries) const
{
//    std::cout << "\nRequested block: (" << b1 << ", " << n1 << "; "
//              << b2 << ", " << n2 << ")" << std::endl;
    if (countAccessedEntries)
        m_accessedEntryCount += n1 * n2;

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
    if (cluster1 && cluster2) {
        // both getdiam2() and dist2() are effectively const, but not declared so
        // in AHMED
        const double diam1 = sqrt(cluster1->getdiam2());
        const double diam2 = sqrt(cluster2->getdiam2());
        const double dist = sqrt(cluster1->extDist2(cluster2));
        minDist = std::max(0., dist - (diam1 + diam2) / 2.);
    }
    // else
        // std::cout << "Warning: clusters not available" << std::endl;
    // std::cout << "minDist: " << minDist << std::endl;

    // This is a non-op for real types. For complex types, it converts a pointer
    // to Ahmed's scomp (resp. dcomp) to a pointer to std::complex<float>
    // (resp. std::complex<double>), which should be perfectly safe since these
    // types have the same binary representation.
    ResultType* data = reinterpret_cast<ResultType*>(ahmedData);

    // Convert AHMED matrix indices into point and DOF indices
    shared_ptr<const ComponentLists> componentLists =
            m_componentListsCache->get(b1, n1);
    shared_ptr<const LocalDofLists<BasisFunctionType> > trialDofLists =
            m_trialDofListsCache->get(b2, n2);

    // Necessary points
    const std::vector<int>& pointIndices =
            componentLists->pointIndices;
    // Necessary components at each point
    const std::vector<std::vector<int> >& componentIndices =
            componentLists->componentIndices;

    typedef typename LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
    // Necessary elements
    const std::vector<int>& trialElementIndices =
            trialDofLists->elementIndices;
    // Necessary local dof indices in each element
    const std::vector<std::vector<LocalDofIndex> >& trialLocalDofs =
            trialDofLists->localDofIndices;
    // Weights of local dofs in each element
    const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights =
            trialDofLists->localDofWeights;
    for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
        for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
            assert(std::abs(trialLocalDofWeights[i][j]) > 0.);
    // Corresponding row and column indices in the matrix to be calculated
    // and stored in ahmedData
    const std::vector<std::vector<int> >& blockRows =
            componentLists->arrayIndices;
    const std::vector<std::vector<int> >& blockCols =
            trialDofLists->arrayIndices;

    arma::Mat<ResultType> result(data, n1, n2, false /*copy_aux_mem*/,
                                 true /*strict*/);
    result.fill(0.);

    if (n2 == 1) {
        // Only one column of the block needed. This means that we need only
        // one local DOF from just one or a few trialElements. Evaluate the
        // local potential operator for one local trial DOF at a time.

        // indices: vector: point index; matrix: component, dof
        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t nTrialElem = 0;
             nTrialElem < trialElementIndices.size();
             ++nTrialElem) {
            const int activeTrialElementIndex = trialElementIndices[nTrialElem];

            // The body of this loop will very probably only run once (single
            // local DOF per trial element)
            for (size_t nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof) {
                LocalDofIndex activeTrialLocalDof =
                        trialLocalDofs[nTrialElem][nTrialDof];
                BasisFunctionType activeTrialLocalDofWeight =
                        trialLocalDofWeights[nTrialElem][nTrialDof];
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
                    m_assemblers[nTerm]->evaluateLocalContributions(
                                pointIndices,
                                activeTrialElementIndex, activeTrialLocalDof,
                                localResult, minDist);
                    for (size_t nPoint = 0;
                         nPoint < pointIndices.size();
                         ++nPoint)
                        for (size_t nComponent = 0;
                             nComponent < componentIndices[nPoint].size();
                             ++nComponent)
                            result(blockRows[nPoint][nComponent], 0) +=
                                    m_termMultipliers[nTerm] *
                                    activeTrialLocalDofWeight *
                                    localResult[nPoint](
                                        componentIndices[nPoint][nComponent], 0);
                }
            }
        }
    }
    else if (n1 == 1) {
        // Only one row of the block needed. This means that we need to
        // evaluate a single component of the local potential operator at
        // a single point.
        assert(pointIndices.size() == 1);
        assert(componentIndices.size() == 1);
        assert(componentIndices[0].size() == 1);

        // indices: vector: trial element; matrix: component, dof
        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
            m_assemblers[nTerm]->evaluateLocalContributions(
                        pointIndices[0], componentIndices[0][0],
                        trialElementIndices, localResult, minDist);
            for (size_t nTrialElem = 0;
                 nTrialElem < trialElementIndices.size();
                 ++nTrialElem)
                for (size_t nTrialDof = 0;
                     nTrialDof < trialLocalDofs[nTrialElem].size();
                     ++nTrialDof)
                    result(0, blockCols[nTrialElem][nTrialDof]) +=
                            m_termMultipliers[nTerm] *
                            trialLocalDofWeights[nTrialElem][nTrialDof] *
                            localResult[nTrialElem](
                                0, trialLocalDofs[nTrialElem][nTrialDof]);
        }
    }
    else { // a "fat" block
        // The whole block or its submatrix needed. This means that we are
        // likely to need all or almost all local DOFs from most elements.
        // Evaluate the full local weak form for each pair of test and trial
        // elements and then select the entries that we need.

        Fiber::_2dArray<arma::Mat<ResultType> > localResult;
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
            m_assemblers[nTerm]->evaluateLocalContributions(
                        pointIndices, trialElementIndices,
                        localResult, minDist);
            for (size_t nTrialElem = 0;
                 nTrialElem < trialElementIndices.size();
                 ++nTrialElem)
                for (size_t nTrialDof = 0;
                     nTrialDof < trialLocalDofs[nTrialElem].size();
                     ++nTrialDof)
                    for (size_t nPoint = 0;
                         nPoint < pointIndices.size();
                         ++nPoint)
                        for (size_t nComponent = 0;
                             nComponent < componentIndices[nPoint].size();
                             ++nComponent)
                            result(blockRows[nPoint][nComponent],
                                   blockCols[nTrialElem][nTrialDof]) +=
                                    m_termMultipliers[nTerm] *
                                    trialLocalDofWeights[nTrialElem][nTrialDof] *
                                    localResult(nPoint, nTrialElem)
                                    (componentIndices[nPoint][nComponent],
                                     trialLocalDofs[nTrialElem][nTrialDof]);
        }
    }
}

template <typename BasisFunctionType, typename ResultType>
void PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::cmpblsym(
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
typename PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::scale(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        const cluster* c1, const cluster* c2) const
{
    //    return m_options.acaOptions().scaling;
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
typename PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::relativeScale(
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
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::
accessedEntryCount() const
{
    return m_accessedEntryCount;
}

template <typename BasisFunctionType, typename ResultType>
void
PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>::
resetAccessedEntryCount()
{
    m_accessedEntryCount = 0;
}

// Explicit instantiations

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(PotentialOperatorAcaAssemblyHelper);

}

#endif // WITH_AHMED
