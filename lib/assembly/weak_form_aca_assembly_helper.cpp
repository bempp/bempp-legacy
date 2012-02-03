#include "weak_form_aca_assembly_helper.hpp"

#include "elementary_linear_operator.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/reverse_index_set.hpp"
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
        const arma::Col<unsigned int>& p2oTestDofs,
        const arma::Col<unsigned int>& p2oTrialDofs,
        const QuadratureSelector<ValueType>& quadSelector,
        const AssemblyOptions& options) :
    m_operator(op), m_view(view),
    m_testSpace(testSpace), m_trialSpace(trialSpace),
    m_p2oTestDofs(p2oTestDofs), m_p2oTrialDofs(p2oTrialDofs),
    m_quadSelector(quadSelector), m_options(options)
{}

template <typename ValueType>
void WeakFormAcaAssemblyHelper<ValueType>::cmpbl(
        unsigned b1, unsigned n1, unsigned b2, unsigned n2, ValueType* data) const
{
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
    findLocalDofs(n1, b1, m_p2oTestDofs, m_testSpace,
                  testElements, testLocalDofs, blockRows);
    findLocalDofs(n2, b2, m_p2oTrialDofs, m_trialSpace,
                  trialElements, trialLocalDofs, blockCols);

    if (n2 == 1)
    {
        // Only one column of the block needed. This means that we need only
        // one local DOF from just one or a few trialElements. Evaluate the
        // local weak form for one local trial DOF at a time.
        arma::Col<ValueType> result(data, n2, false /*copy_aux_mem*/,
                                    true /*strict*/);
        result.zeros();

        std::vector<const EntityPointer<0>*> activeTrialElements(1);
        const EntityPointer<0>*& activeTrialElement = activeTrialElements[0];
        std::vector<arma::Col<ValueType> > localResult;
        for (int nTrialElem = 0;
             nTrialElem < trialElements.size();
             ++nTrialElem)
        {
            activeTrialElement = trialElements[nTrialElem];
            for (int nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof)
            {
                LocalDofIndex activeTrialLocalDof =
                        trialLocalDofs[nTrialElem][nTrialDof];
                m_operator.evaluateLocalWeakForms(
                            testElements,
                            activeTrialElement, activeTrialLocalDof,
                            m_testSpace, m_trialSpace, m_quadSelector,
                            m_options, localResult);
                for (int nTestElem = 0;
                     nTestElem < testElements.size();
                     ++nTestElem)
                {
                    for (int nTestDof = 0;
                         nTestDof < testLocalDofs[nTestElem].size();
                         ++nTestDof)
                        result(blockRows[nTestElem][nTestDof]) +=
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

        // TODO: code simlar to the above, with trial and test elements swapped.
    }
    else // a "fat" block
    {
        // The whole block or its submatrix needed. This means that we are
        // likely to need all or almost all local DOFs from most elements.
        // Evaluate the full local weak form for each pair of test and trial
        // elements and then select the entries that we need.
        arma::Mat<ValueType> result(data, n1, n2, false /*copy_aux_mem*/,
                                         true /*strict*/);
        result.zeros();

        Array2D<arma::Mat<ValueType> > localResult;
        m_operator.evaluateLocalWeakForms(
                    testElements, trialElements,
                    m_testSpace, m_trialSpace, m_quadSelector, m_options,
                    localResult);
        for (int nTrialElem = 0;
             nTrialElem < trialElements.size();
             ++nTrialElem)
        {
            for (int nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size();
                 ++nTrialDof)
            {
                for (int nTestElem = 0;
                     nTestElem < testElements.size();
                     ++nTestElem)
                {
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
        }
    }
}

template <typename ValueType>
void WeakFormAcaAssemblyHelper<ValueType>::findLocalDofs(
        int start,
        int length,
        const arma::Col<unsigned int>& p2o,
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

    vector<GlobalDofIndex> dofs(length);
    for (int i = 0; i < length; ++i)
        dofs[i] = p2o(start + i);

    vector<vector<LocalDof> > localDofs;
    space.global2localDofs(dofs, localDofs);

    typedef set<pair<LocalDofIndex, int> > LocalDofSet;
    typedef map<EntityIndex, LocalDofSet> LocalDofMap;
    LocalDofMap requiredLocalDofs;
    // arrayIndex: row or column in the matrix "result" to be filled
    for (int arrayIndex = 0; arrayIndex < localDofs.size(); ++arrayIndex)
    {
        const vector<LocalDof>& currentLocalDofs = localDofs[arrayIndex];
        for (int j = 0; j < currentLocalDofs.size(); ++j)
            requiredLocalDofs[currentLocalDofs[j].entityIndex]
                    .insert(make_pair(currentLocalDofs[j].dofIndex, arrayIndex));
    }

    const int elementCount = requiredLocalDofs.size();
    // vector<EntityIndex> elementIndices;
    // elementIndices.reserve(elementCount);
    elements.reserve(elementCount);
    localDofIndices.reserve(elementCount);
    arrayIndices.reserve(elementCount);

    const ReverseIndexSet& reverseIndexSet = m_view.reverseIndexSet();

    int e = 0;
    for (LocalDofMap::const_iterator mapIt = requiredLocalDofs.begin();
         mapIt != requiredLocalDofs.end(); ++mapIt, ++e)
    {
        // elementIndices[e] = mapIt->first;
        elements[e] = &reverseIndexSet.entityPointerCodim0(mapIt->first);
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
