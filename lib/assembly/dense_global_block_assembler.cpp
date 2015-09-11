// Copyright (C) 2011-2015 by the BEM++ Authors
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

#include "dense_global_block_assembler.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "assembly_options.hpp"
#include "evaluation_options.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "context.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/reverse_element_mapper.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/complex_aux.hpp"
#include "../common/eigen_support.hpp"
#include <stdexcept>
#include <iostream>

#include <tbb/parallel_for_each.h>
#include <tbb/spin_mutex.h>
//#include <tbb/tick_count.h>

#include <unordered_map>
#include <utility>

namespace Bempp {

// Helper functions and classes

namespace {

// Body of parallel loop

template <typename BasisFunctionType, typename ResultType>
class DenseWeakFormAssemblerLoopBody {
public:
  typedef tbb::spin_mutex MutexType;

  DenseWeakFormAssemblerLoopBody(
      int rowStart, int colStart, const std::vector<int> &testIndices,
      const std::unordered_map<int, std::vector<GlobalDofIndex>>
          &testGlobalDofs,
      const std::unordered_map<int, std::vector<GlobalDofIndex>>
          &trialGlobalDofs,
      const std::unordered_map<int, std::vector<BasisFunctionType>>
          &testLocalDofWeights,
      const std::unordered_map<int, std::vector<BasisFunctionType>>
          &trialLocalDofWeights,
      Fiber::LocalAssemblerForIntegralOperators<ResultType> &assembler,
      Matrix<ResultType> &result, MutexType &mutex)
      : m_rowStart(rowStart), m_colStart(colStart),
        m_testIndices(testIndices),
        m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs),
        m_testLocalDofWeights(testLocalDofWeights),
        m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler),
        m_result(result), m_mutex(mutex) {}

  void
  operator()(const std::pair<const int, std::vector<GlobalDofIndex>> &p) const {
    std::vector<Matrix<ResultType>> localResult;
    int trialElementIndex = p.first;
    const std::vector<GlobalDofIndex> &trialElementDofs = p.second;
    const std::vector<BasisFunctionType> &trialElementWeights =
        m_trialLocalDofWeights.at(trialElementIndex);
    // Evaluate integrals over pairs of the current trial element and
    // all the test elements
    m_assembler.evaluateLocalWeakForms(
        TEST_TRIAL, m_testIndices, trialElementIndex, ALL_DOFS, localResult);

    {
      MutexType::scoped_lock lock(m_mutex);
      int count = 0;
      for (const auto &testPair : m_testGlobalDofs) {
        int testElementIndex = testPair.first;
        const std::vector<GlobalDofIndex> &testElementDofs = testPair.second;
        const std::vector<BasisFunctionType> &testElementWeights =
            m_testLocalDofWeights.at(testElementIndex);
        for (int trialDof = 0; trialDof < trialElementDofs.size(); ++trialDof) {
          int trialGlobalDof = trialElementDofs[trialDof];
          if (trialGlobalDof < 0)
            continue;

          for (int testDof = 0; testDof < testElementDofs.size(); ++testDof) {
            int testGlobalDof = testElementDofs[testDof];
            if (testGlobalDof < 0)
              continue;

            m_result(testGlobalDof - m_rowStart, trialGlobalDof - m_colStart) +=
                conj(testElementWeights[testDof]) *
                trialElementWeights[trialDof] *
                localResult[count](testDof, trialDof);
          }
        }
        count++;
      }
    }
  }

private:
  int m_rowStart, m_colStart;
  const std::vector<int> &m_testIndices;
  const std::unordered_map<int, std::vector<GlobalDofIndex>> &m_testGlobalDofs;
  const std::unordered_map<int, std::vector<GlobalDofIndex>> &m_trialGlobalDofs;
  const std::unordered_map<int, std::vector<BasisFunctionType>>
      &m_testLocalDofWeights;
  const std::unordered_map<int, std::vector<BasisFunctionType>>
      &m_trialLocalDofWeights;
  // mutable OK because Assembler is thread-safe. (Alternative to "mutable"
  // here:
  // make assembler's internal integrator map mutable)
  typename Fiber::LocalAssemblerForIntegralOperators<ResultType> &m_assembler;
  // mutable OK because write access to this matrix is protected by a mutex
  Matrix<ResultType> &m_result;

  // mutex must be mutable because we need to lock and unlock it
  MutexType &m_mutex;
};

template <typename BasisFunctionType>
void gatherElementInformation(
    int start, int end, const Space<BasisFunctionType> &space,
    std::unordered_map<int, std::vector<GlobalDofIndex>> &indexMap,
    std::unordered_map<int, std::vector<BasisFunctionType>> &dofWeights) {
  auto gridView = space.grid()->leafView();
  std::vector<GlobalDofIndex> globalDofs;
  for (int i = start; i < end; ++i)
    globalDofs.push_back(i);

  std::vector<std::vector<LocalDof>> global2localDofs;
  std::vector<std::vector<BasisFunctionType>> weights; // Not needed
  space.global2localDofs(globalDofs, global2localDofs, weights);

  const ReverseElementMapper& reverseMapper = gridView->reverseElementMapper();

  for (const auto &localDofs : global2localDofs)
    for (const auto &localDof : localDofs) {
      int elementIndex = localDof.entityIndex;
      if (!indexMap.count(elementIndex)) {
        const auto& element = reverseMapper.entityPointer(elementIndex).entity();
        indexMap[elementIndex] = std::vector<GlobalDofIndex>();
        dofWeights[elementIndex] = std::vector<BasisFunctionType>();
        space.getGlobalDofs(element, indexMap.at(elementIndex),
                            dofWeights.at(elementIndex));
        for (auto &dof : indexMap.at(elementIndex))
          if (dof < start || dof >= end)
            dof = -1;
      }
    }
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
assembleDenseBlock(
    int rowStart, int rowEnd, int colStart, int colEnd,
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler,
    const ParameterList &parameterList) {

  int numberOfRows = rowEnd - rowStart;
  int numberOfColumns = colEnd - colStart;

  if (colEnd > trialSpace.globalDofCount() ||
      rowEnd > testSpace.globalDofCount() || colStart < 0 || rowStart < 0)
    throw std::runtime_error("DenseGlobalBlockAssember::assembleWeakForm(): "
                             "Indices out of bounds");

  Context<BasisFunctionType, ResultType> context(parameterList);
  const AssemblyOptions &options = context.assemblyOptions();

  // Create the operator's matrix
  Matrix<ResultType> result(numberOfRows, numberOfColumns);
  result.setZero();

  std::unordered_map<int, std::vector<GlobalDofIndex>> trialIndexMap,
      testIndexMap;
  std::unordered_map<int, std::vector<BasisFunctionType>> trialDofWeights,
      testDofWeights;

  gatherElementInformation(colStart, colEnd, trialSpace, trialIndexMap,
                           trialDofWeights);
  gatherElementInformation(rowStart, rowEnd, testSpace, testIndexMap,
                           testDofWeights);

  std::vector<int> testIndices;
  testIndices.reserve(testIndexMap.size());
  for (const auto &p : testIndexMap)
    testIndices.push_back(p.first);

  typedef DenseWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
  typename Body::MutexType mutex;

  {
    Fiber::SerialBlasRegion region;
    tbb::parallel_for_each(trialIndexMap.begin(), trialIndexMap.end(),
                           Body(rowStart, colStart, testIndices, testIndexMap,
                                trialIndexMap, testDofWeights, trialDofWeights,
                                assembler, result, mutex));
  }

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return shared_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT>>  \
    assembleDenseBlock(int, int, int, int, \
            const Space<BASIS>&, const Space<BASIS>&, \
            Fiber::LocalAssemblerForIntegralOperators<RESULT>&, \
            const ParameterList&)

FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FREE_FUNCTIONS);

} // namespace Bempp
