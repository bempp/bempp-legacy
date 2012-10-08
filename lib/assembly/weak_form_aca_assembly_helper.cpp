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

#include "discrete_boundary_operator.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../fiber/types.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/reverse_element_mapper.hpp"
#include "../space/space.hpp"
#include "assembly_options.hpp"

#include "ahmed_complex.hpp"

#include <map>
#include <set>
#include <utility>
#include <cstdio>

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
    m_indexWithGlobalDofs(m_options.acaOptions().globalAssemblyBeforeCompression)
{
    if (!m_indexWithGlobalDofs && !m_sparseTermsToAdd.empty())
        throw std::invalid_argument(
                "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
                "combining sparse and dense terms with "
                "globalAssemblyBeforeCompression set to false "
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
}

template <typename BasisFunctionType, typename ResultType>
void WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::cmpbl(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2,
        AhmedResultType* ahmedData) const
{
//    std::cout << "\nRequested block: (" << b1 << ", " << n1 << "; "
//              << b2 << ", " << n2 << ")" << std::endl;

    // This is a non-op for real types. For complex types, it converts a pointer
    // to Ahmed's scomp (resp. dcomp) to a pointer to std::complex<float>
    // (resp. std::complex<double>), which should be perfectly safe since these
    // types have the same binary representation.
    ResultType* data = reinterpret_cast<ResultType*>(ahmedData);

    // Requested original matrix indices
    std::vector<DofIndex> testOriginalIndices;
    std::vector<DofIndex> trialOriginalIndices;
    // Necessary elements
    std::vector<int> testElementIndices;
    std::vector<int> trialElementIndices;
    // Necessary local dof indices in each element
    std::vector<std::vector<LocalDofIndex> > testLocalDofs;
    std::vector<std::vector<LocalDofIndex> > trialLocalDofs;
    // Corresponding row and column indices in the matrix to be calculated
    // and stored in ahmedData
    std::vector<std::vector<int> > blockRows;
    std::vector<std::vector<int> > blockCols;

    // Fill the above arrays
    findLocalDofs(b1, n1, m_p2oTestDofs, m_testSpace,
                  testOriginalIndices, testElementIndices, testLocalDofs,
                  blockRows);
    findLocalDofs(b2, n2, m_p2oTrialDofs, m_trialSpace,
                  trialOriginalIndices, trialElementIndices, trialLocalDofs,
                  blockCols);

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
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
                {
                    m_assemblers[nTerm]->evaluateLocalWeakForms(
                                Fiber::TEST_TRIAL, testElementIndices,
                                activeTrialElementIndex, activeTrialLocalDof,
                                localResult);
                    for (size_t nTestElem = 0;
                         nTestElem < testElementIndices.size();
                         ++nTestElem)
                        for (size_t nTestDof = 0;
                             nTestDof < testLocalDofs[nTestElem].size();
                             ++nTestDof)
                            result(blockRows[nTestElem][nTestDof], 0) +=
                                    m_denseTermsMultipliers[nTerm] *
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
                for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
                {
                    m_assemblers[nTerm]->evaluateLocalWeakForms(
                                Fiber::TRIAL_TEST, trialElementIndices,
                                activeTestElementIndex, activeTestLocalDof,
                                localResult);
                    for (size_t nTrialElem = 0;
                         nTrialElem < trialElementIndices.size();
                         ++nTrialElem)
                        for (size_t nTrialDof = 0;
                             nTrialDof < trialLocalDofs[nTrialElem].size();
                             ++nTrialDof)
                            result(0, blockCols[nTrialElem][nTrialDof]) +=
                                    m_denseTermsMultipliers[nTerm] *
                                    localResult[nTrialElem](trialLocalDofs[nTrialElem][nTrialDof]);
                }
            }
        }
    }
    else // a "fat" block
    {
        // The whole block or its submatrix needed. This means that we are
        // likely to need all or almost all local DOFs from most elements.
        // Evaluate the full local weak form for each pair of test and trial
        // elements and then select the entries that we need.

        Fiber::_2dArray<arma::Mat<ResultType> > localResult;
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
        {
            m_assemblers[nTerm]->evaluateLocalWeakForms(
                        testElementIndices, trialElementIndices, localResult);
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
                                    localResult(nTestElem, nTrialElem)
                                    (testLocalDofs[nTestElem][nTestDof],
                                     trialLocalDofs[nTrialElem][nTrialDof]);
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
        unsigned b1, unsigned n1, AhmedResultType* ahmedData) const
{
    // Calculate results as usual (without taking symmetry into account)
    arma::Mat<ResultType> block(n1, n1);
    cmpbl(b1, n1, b1, n1, ahmedCast(block.memptr()));

    // and now store the upper part of the matrix in the memory block
    // provided by Ahmed
    for (size_t col = 0; col < n1; ++col)
        for (size_t row = 0; row <= col; ++row)
            *ahmedData++ = ahmedCast(block(row, col));
}

template <typename BasisFunctionType, typename ResultType>
typename WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::MagnitudeType
WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::scale(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2) const
{
    return m_options.acaOptions().scaling;
}

template <typename BasisFunctionType, typename ResultType>
void WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>::
findLocalDofs(
        int start,
        int indexCount,
        const std::vector<unsigned int>& p2o,
        const Space<BasisFunctionType>& space,
        std::vector<DofIndex>& originalIndices,
        std::vector<int>& elementIndices,
        std::vector<std::vector<LocalDofIndex> >& localDofIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    using std::make_pair;
    using std::map;
    using std::pair;
    using std::set;
    using std::vector;

    // Convert permuted indices into original indices
    originalIndices.resize(indexCount);
    for (int i = 0; i < indexCount; ++i)
        originalIndices[i] = p2o[start + i];

    // set of pairs (local dof index, array index)
    typedef set<pair<LocalDofIndex, int> > LocalDofSet;
    typedef map<EntityIndex, LocalDofSet> LocalDofMap;

    // Temporary map: entityIndex -> set(localDofIndex, arrayIndex)
    // with arrayIndex standing for the index of the row or column in the matrix
    // that needs to be returned to Ahmed.
    LocalDofMap requiredLocalDofs;

    // Retrieve lists of local DOFs corresponding to original indices,
    // treated either as global DOFs (if m_indexWithGlobalDofs is true)
    // or flat local DOFs (if m_indexWithGlobalDofs is false)
    if (m_indexWithGlobalDofs)
    {
        vector<vector<LocalDof> > localDofs;
        space.global2localDofs(originalIndices, localDofs);

        for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex)
        {
            const vector<LocalDof>& currentLocalDofs = localDofs[arrayIndex];
            for (size_t j = 0; j < currentLocalDofs.size(); ++j)
                requiredLocalDofs[currentLocalDofs[j].entityIndex]
                        .insert(make_pair(currentLocalDofs[j].dofIndex, arrayIndex));
        }
    }
    else
    {
        vector<LocalDof> localDofs;
        space.flatLocal2localDofs(originalIndices, localDofs);

        for (int arrayIndex = 0; arrayIndex < indexCount; ++arrayIndex)
        {
            const LocalDof& currentLocalDof = localDofs[arrayIndex];
            requiredLocalDofs[currentLocalDof.entityIndex]
                    .insert(make_pair(currentLocalDof.dofIndex, arrayIndex));
        }
    }

    // Use the temporary map requiredLocalDofs to build the three output vectors
    const int elementCount = requiredLocalDofs.size();
    // vector<EntityIndex> elementIndices;
    // elementIndices.reserve(elementCount);

    elementIndices.resize(elementCount);
    localDofIndices.clear();
    localDofIndices.resize(elementCount);
    arrayIndices.clear();
    arrayIndices.resize(elementCount);

    // const ReverseElementMapper& mapper = m_view.reverseElementMapper();

    int e = 0;
    for (LocalDofMap::const_iterator mapIt = requiredLocalDofs.begin();
         mapIt != requiredLocalDofs.end(); ++mapIt, ++e)
    {
        elementIndices[e] = mapIt->first;
//        elements[e] = &mapper.entityPointer(mapIt->first);
        for (LocalDofSet::const_iterator setIt = mapIt->second.begin();
             setIt != mapIt->second.end(); ++setIt)
        {
            localDofIndices[e].push_back(setIt->first);
            arrayIndices[e].push_back(setIt->second);
        }
    }
}

// Explicit instantiations

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(WeakFormAcaAssemblyHelper);

}

#endif // WITH_AHMED
