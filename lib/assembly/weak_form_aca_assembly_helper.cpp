#include "weak_form_aca_assembly_helper.hpp"

#include "elementary_linear_operator.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../fiber/types.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/reverse_element_mapper.hpp"
#include "../space/space.hpp"

#include <map>
#include <set>
#include <utility>

namespace Bempp
{

template <typename ValueType>
WeakFormAcaAssemblyHelper<ValueType>::WeakFormAcaAssemblyHelper(
        const ElementaryLinearOperator<ValueType>& op,
        const GridView& view,
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const std::vector<unsigned int>& p2oTestDofs,
        const std::vector<unsigned int>& p2oTrialDofs,
        IntegrationManager& intMgr,
        const AssemblyOptions& options) :
    m_operator(op), m_view(view),
    m_testSpace(testSpace), m_trialSpace(trialSpace),
    m_p2oTestDofs(p2oTestDofs), m_p2oTrialDofs(p2oTrialDofs),
    m_intMgr(intMgr), m_options(options)
{}

template <typename ValueType>
void WeakFormAcaAssemblyHelper<ValueType>::cmpbl(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2, ValueType* data) const
{
//    std::cout << "\nRequested block: (" << b1 << ", " << n1 << "; "
//              << b2 << ", " << n2 << ")" << std::endl;

    // Necessary elements
    std::vector<const EntityPointer<0>*> testElements;
    std::vector<const EntityPointer<0>*> trialElements;
    // Necessary local dof indices in each element
    std::vector<std::vector<LocalDofIndex> > testLocalDofs;
    std::vector<std::vector<LocalDofIndex> > trialLocalDofs;
    // Corresponding row and column indices in the matrix to be calculated
    std::vector<std::vector<int> > blockRows;
    std::vector<std::vector<int> > blockCols;

    // Fill the above arrays
    findLocalDofs(b1, n1, m_p2oTestDofs, m_testSpace,
                  testElements, testLocalDofs, blockRows);
    findLocalDofs(b2, n2, m_p2oTrialDofs, m_trialSpace,
                  trialElements, trialLocalDofs, blockCols);

    if (n2 == 1)
    {
        // Only one column of the block needed. This means that we need only
        // one local DOF from just one or a few trialElements. Evaluate the
        // local weak form for one local trial DOF at a time.
        arma::Col<ValueType> result(data, n1, false /*copy_aux_mem*/,
                                    true /*strict*/);
        result.fill(0.);

        std::vector<arma::Mat<ValueType> > localResult;
        for (int nTrialElem = 0;
             nTrialElem < trialElements.size();
             ++nTrialElem)
        {
            const EntityPointer<0>& activeTrialElement = *trialElements[nTrialElem];
            // The body of this loop will very probably only run once (single
            // local DOF per trial element)
            for (int nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof)
            {
                LocalDofIndex activeTrialLocalDof =
                        trialLocalDofs[nTrialElem][nTrialDof];
                m_operator.evaluateLocalWeakForms(
                            Fiber::TEST_TRIAL,
                            testElements, activeTrialElement, activeTrialLocalDof,
                            m_testSpace, m_trialSpace, m_intMgr, localResult);
                for (int nTestElem = 0;
                     nTestElem < testElements.size();
                     ++nTestElem)
                    for (int nTestDof = 0;
                         nTestDof < testLocalDofs[nTestElem].size();
                         ++nTestDof)
                        result(blockRows[nTestElem][nTestDof]) +=
                                localResult[nTestElem](testLocalDofs[nTestElem][nTestDof]);
            }
        }
    }
    else if (n1 == 1) // very few testElements
    {
        // Only one row of the block needed. This means that we need only
        // one local DOF from just one or a few testElements. Evaluate the
        // local weak form for one local test DOF at a time.

        arma::Row<ValueType> result(data, n2, false /*copy_aux_mem*/,
                                    true /*strict*/);
        result.fill(0.);

        std::vector<arma::Mat<ValueType> > localResult;
        for (int nTestElem = 0;
             nTestElem < testElements.size();
             ++nTestElem)
        {
            const EntityPointer<0>& activeTestElement = *testElements[nTestElem];
            // The body of this loop will very probably only run once (single
            // local DOF per test element)
            for (int nTestDof = 0;
                 nTestDof < testLocalDofs[nTestElem].size();
                 ++nTestDof)
            {
                LocalDofIndex activeTestLocalDof =
                        testLocalDofs[nTestElem][nTestDof];
                m_operator.evaluateLocalWeakForms(
                            Fiber::TRIAL_TEST,
                            trialElements, activeTestElement, activeTestLocalDof,
                            m_trialSpace, m_testSpace, m_intMgr, localResult);
                for (int nTrialElem = 0;
                     nTrialElem < trialElements.size();
                     ++nTrialElem)
                    for (int nTrialDof = 0;
                         nTrialDof < trialLocalDofs[nTrialElem].size();
                         ++nTrialDof)
                        result(blockCols[nTrialElem][nTrialDof]) +=
                                localResult[nTrialElem](trialLocalDofs[nTrialElem][nTrialDof]);
            }
        }
    }
    else // a "fat" block
    {
        // The whole block or its submatrix needed. This means that we are
        // likely to need all or almost all local DOFs from most elements.
        // Evaluate the full local weak form for each pair of test and trial
        // elements and then select the entries that we need.
        arma::Mat<ValueType> result(data, n1, n2, false /*copy_aux_mem*/,
                                    true /*strict*/);
        result.fill(0.);

        Fiber::Array2D<arma::Mat<ValueType> > localResult;
        m_operator.evaluateLocalWeakForms(
                    testElements, trialElements,
                    m_testSpace, m_trialSpace, m_intMgr, localResult);
        for (int nTrialElem = 0;
             nTrialElem < trialElements.size();
             ++nTrialElem)
            for (int nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof)
                for (int nTestElem = 0;
                     nTestElem < testElements.size();
                     ++nTestElem)
                    for (int nTestDof = 0;
                         nTestDof < testLocalDofs[nTestElem].size();
                         ++nTestDof)
                        result(blockRows[nTestElem][nTestDof],
                               blockCols[nTrialElem][nTrialDof]) +=
                                localResult(nTestElem, nTrialElem)
                                (testLocalDofs[nTestElem][nTestDof],
                                 trialLocalDofs[nTrialElem][nTrialDof]);
    }
}

template <typename ValueType>
ValueType WeakFormAcaAssemblyHelper<ValueType>::scale(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2) const
{
    return 1.;
}

template <typename ValueType>
void WeakFormAcaAssemblyHelper<ValueType>::findLocalDofs(
        int start,
        int globalDofCount,
        const std::vector<unsigned int>& p2o,
        const Space<ValueType>& space,
        std::vector<const EntityPointer<0>*>& elements,
        std::vector<std::vector<LocalDofIndex> >& localDofIndices,
        std::vector<std::vector<int> >& arrayIndices) const
{
    using std::make_pair;
    using std::map;
    using std::pair;
    using std::set;
    using std::vector;

    // Convert permuted indices into the original global DOF indices
    vector<GlobalDofIndex> dofs(globalDofCount);
    for (int i = 0; i < globalDofCount; ++i)
        dofs[i] = p2o[start + i];

    // Retrieve lists of local DOFs corresponding to the global DOFs
    vector<vector<LocalDof> > localDofs;
    space.global2localDofs(dofs, localDofs);

    // set of pairs (local dof index, array index)
    typedef set<pair<LocalDofIndex, int> > LocalDofSet;
    typedef map<EntityIndex, LocalDofSet> LocalDofMap;
    // Temporary map: entityIndex -> set(localDofIndex, arrayIndex)
    // with arrayIndex standing for the index of the row or column in the matrix
    // that needs to be returned to Ahmed.
    LocalDofMap requiredLocalDofs;
    for (int arrayIndex = 0; arrayIndex < globalDofCount; ++arrayIndex)
    {
        const vector<LocalDof>& currentLocalDofs = localDofs[arrayIndex];
        for (int j = 0; j < currentLocalDofs.size(); ++j)
            requiredLocalDofs[currentLocalDofs[j].entityIndex]
                    .insert(make_pair(currentLocalDofs[j].dofIndex, arrayIndex));
    }

    // Use the temporary map requiredLocalDofs to build the three output vectors
    const int elementCount = requiredLocalDofs.size();
    // vector<EntityIndex> elementIndices;
    // elementIndices.reserve(elementCount);

    elements.resize(elementCount);
    localDofIndices.clear();
    localDofIndices.resize(elementCount);
    arrayIndices.clear();
    arrayIndices.resize(elementCount);

    const ReverseElementMapper& mapper = m_view.reverseElementMapper();

    int e = 0;
    for (LocalDofMap::const_iterator mapIt = requiredLocalDofs.begin();
         mapIt != requiredLocalDofs.end(); ++mapIt, ++e)
    {
        // elementIndices[e] = mapIt->first;
        elements[e] = &mapper.entityPointer(mapIt->first);
        for (LocalDofSet::const_iterator setIt = mapIt->second.begin();
             setIt != mapIt->second.end(); ++setIt)
        {
            localDofIndices[e].push_back(setIt->first);
            arrayIndices[e].push_back(setIt->second);
        }
    }
}

// Explicit instantiations

#ifdef COMPILE_FOR_FLOAT
template class WeakFormAcaAssemblyHelper<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class WeakFormAcaAssemblyHelper<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class WeakFormAcaAssemblyHelper<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class WeakFormAcaAssemblyHelper<std::complex<double> >;
#endif

}
