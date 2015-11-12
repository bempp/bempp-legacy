
#include "potential_operator_hmat_assembly_helper.hpp"

#include "../space/space.hpp"
#include "component_lists_cache.hpp"
#include "local_dof_lists_cache.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
PotentialOperatorHMatAssemblyHelper<BasisFunctionType, ResultType>::
    PotentialOperatorHMatAssemblyHelper(
        const Matrix<CoordinateType> &points,
        const Space<BasisFunctionType> &trialSpace,
        const shared_ptr<const hmat::DefaultBlockClusterTreeType>
            &blockClusterTree,
        LocalAssembler &assembler, const ParameterList &parameterList)
    : m_points(points), m_trialSpace(trialSpace),
      m_blockClusterTree(blockClusterTree), m_assembler(assembler),
      m_trialDofListsCache(new LocalDofListsCache<BasisFunctionType>(
          m_trialSpace,
          blockClusterTree->columnClusterTree()->hMatDofToOriginalDofMap(),
          true)),
      m_parameterList(parameterList),
      m_componentCount(assembler.resultDimension()),
      m_componentListsCache(new ComponentListsCache(
          blockClusterTree->rowClusterTree()->hMatDofToOriginalDofMap(),
          m_componentCount)) {}

template <typename BasisFunctionType, typename ResultType>
typename PotentialOperatorHMatAssemblyHelper<BasisFunctionType,
                                             ResultType>::MagnitudeType
PotentialOperatorHMatAssemblyHelper<BasisFunctionType, ResultType>::
    estimateMinimumDistance(const hmat::DefaultBlockClusterTreeNodeType
                                &blockClusterTreeNode) const {

  MagnitudeType dist =
      MagnitudeType(blockClusterTreeNode.data()
                        .rowClusterTreeNode->data()
                        .boundingBox.distance(blockClusterTreeNode.data()
                                                  .columnClusterTreeNode->data()
                                                  .boundingBox));
  return dist;
}

template <typename BasisFunctionType, typename ResultType>
void PotentialOperatorHMatAssemblyHelper<BasisFunctionType, ResultType>::
    computeMatrixBlock(
        const hmat::IndexRangeType &testIndexRange,
        const hmat::IndexRangeType &trialIndexRange,
        const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
        Matrix<ResultType> &data) const {

  auto numberOfTestIndices = testIndexRange[1] - testIndexRange[0];
  auto numberOfTrialIndices = trialIndexRange[1] - trialIndexRange[0];


  const CoordinateType minDist = estimateMinimumDistance(blockClusterTreeNode);

  shared_ptr<const ComponentLists> componentLists =
      m_componentListsCache->get(testIndexRange[0], numberOfTestIndices);
  shared_ptr<const LocalDofLists<BasisFunctionType>> trialDofLists =
      m_trialDofListsCache->get(trialIndexRange[0], numberOfTrialIndices);

  // Necessary points
  const std::vector<int> &pointIndices = componentLists->pointIndices;
  // Necessary components at each point
  const std::vector<std::vector<int>> &componentIndices =
      componentLists->componentIndices;

  typedef typename LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
  // Necessary elements
  const std::vector<int> &trialElementIndices = trialDofLists->elementIndices;
  // Necessary local dof indices in each element
  const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs =
      trialDofLists->localDofIndices;
  // Weights of local dofs in each element
  const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights =
      trialDofLists->localDofWeights;
  for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
    for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
      assert(std::abs(trialLocalDofWeights[i][j]) > 0.);
  // Corresponding row and column indicies in the matrix to be calculated.
  const std::vector<std::vector<int>> &blockRows = componentLists->arrayIndices;
  const std::vector<std::vector<int>> &blockCols = trialDofLists->arrayIndices;

  data.resize(numberOfTestIndices, numberOfTrialIndices);
  data.setZero();

  if (numberOfTrialIndices == 1) {
    // Only one column of the block needed. This means that we need only
    // one local DOF from just one or a few trialElements. Evaluate the
    // local potential operator for one local trial DOF at a time.

    // indices: vector: point index; matrix: component, dof
    std::vector<Matrix<ResultType>> localResult;
    for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
         ++nTrialElem) {
      const int activeTrialElementIndex = trialElementIndices[nTrialElem];

      // The body of this loop will very probably only run once (single
      // local DOF per trial element)
      for (size_t nTrialDof = 0; nTrialDof < trialLocalDofs[nTrialElem].size();
           ++nTrialDof) {
        LocalDofIndex activeTrialLocalDof =
            trialLocalDofs[nTrialElem][nTrialDof];
        BasisFunctionType activeTrialLocalDofWeight =
            trialLocalDofWeights[nTrialElem][nTrialDof];
        m_assembler.evaluateLocalContributions(
            pointIndices, activeTrialElementIndex, activeTrialLocalDof,
            localResult, minDist);
        for (size_t nPoint = 0; nPoint < pointIndices.size(); ++nPoint)
          for (size_t nComponent = 0;
               nComponent < componentIndices[nPoint].size(); ++nComponent)
            data(blockRows[nPoint][nComponent], 0) +=
                activeTrialLocalDofWeight *
                localResult[nPoint](componentIndices[nPoint][nComponent], 0);
      }
    }
  } else if (numberOfTestIndices == 1) {
    // Only one row of the block needed. This means that we need to
    // evaluate a single component of the local potential operator at
    // a single point.
    assert(pointIndices.size() == 1);
    assert(componentIndices.size() == 1);
    assert(componentIndices[0].size() == 1);

    // indices: vector: trial element; matrix: component, dof
    std::vector<Matrix<ResultType>> localResult;
    m_assembler.evaluateLocalContributions(
        pointIndices[0], componentIndices[0][0], trialElementIndices,
        localResult, minDist);
    for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
         ++nTrialElem)
      for (size_t nTrialDof = 0; nTrialDof < trialLocalDofs[nTrialElem].size();
           ++nTrialDof)
        data(0, blockCols[nTrialElem][nTrialDof]) +=
            trialLocalDofWeights[nTrialElem][nTrialDof] *
            localResult[nTrialElem](0, trialLocalDofs[nTrialElem][nTrialDof]);
  } else { // a "fat" block
    // The whole block or its submatrix needed. This means that we are
    // likely to need all or almost all local DOFs from most elements.
    // Evaluate the full local weak form for each pair of test and trial
    // elements and then select the entries that we need.

    Fiber::_2dArray<Matrix<ResultType>> localResult;
    m_assembler.evaluateLocalContributions(pointIndices, trialElementIndices,
                                            localResult, minDist);
    for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
         ++nTrialElem)
      for (size_t nTrialDof = 0; nTrialDof < trialLocalDofs[nTrialElem].size();
           ++nTrialDof)
        for (size_t nPoint = 0; nPoint < pointIndices.size(); ++nPoint)
          for (size_t nComponent = 0;
               nComponent < componentIndices[nPoint].size(); ++nComponent)
            data(blockRows[nPoint][nComponent],
                   blockCols[nTrialElem][nTrialDof]) +=
                trialLocalDofWeights[nTrialElem][nTrialDof] *
                localResult(nPoint,
                            nTrialElem)(componentIndices[nPoint][nComponent],
                                        trialLocalDofs[nTrialElem][nTrialDof]);
  }
}

template <typename BasisFunctionType, typename ResultType>
double
PotentialOperatorHMatAssemblyHelper<BasisFunctionType, ResultType>::scale(
    const hmat::DefaultBlockClusterTreeNodeType &node) const {

  MagnitudeType dist = this->estimateMinimumDistance(node);
  return static_cast<double>(m_assembler.estimateRelativeScale(dist));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    PotentialOperatorHMatAssemblyHelper);
}
