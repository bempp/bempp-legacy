// Copyright (C) 2011-2013 by the BEM++ Authors
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

#include "dense_global_assembler.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "assembly_options.hpp"
#include "evaluation_options.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "context.hpp"

#include "../common/auto_timer.hpp"
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
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/complex_aux.hpp"
#include "../common/eigen_support.hpp"
#include <stdexcept>
#include <iostream>
#include <chrono>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
//#include <tbb/tick_count.h>

namespace Bempp {

// Helper functions and classes

namespace {

// Body of parallel loop

template <typename BasisFunctionType, typename ResultType>
class DenseWeakFormAssemblerLoopBody {
public:
  typedef tbb::spin_mutex MutexType;

  DenseWeakFormAssemblerLoopBody(
      const std::vector<int> &testIndices,
      const std::vector<std::vector<GlobalDofIndex>> &testGlobalDofs,
      const std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
      const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
      const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
      Fiber::LocalAssemblerForIntegralOperators<ResultType> &assembler,
      Matrix<ResultType> &result, Matrix<MutexType> &mutex)
      : m_testIndices(testIndices), m_testGlobalDofs(testGlobalDofs),
        m_trialGlobalDofs(trialGlobalDofs),
        m_testLocalDofWeights(testLocalDofWeights),
        m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler),
        m_result(result), m_mutex(mutex) {}

  void operator()(const tbb::blocked_range<int> &r) const {
    const int testElementCount = m_testIndices.size();
    std::vector<Matrix<ResultType>> localResult;
    for (int trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
      // Handle this trial element only if it contributes to any global DOFs.
      bool skipTrialElement = true;
      const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
      for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
        int trialGlobalDof = m_trialGlobalDofs[trialIndex][trialDof];
        if (trialGlobalDof >= 0) {
          skipTrialElement = false;
          break;
        }
      }
      if (skipTrialElement)
        continue;

      // Evaluate integrals over pairs of the current trial element and
      // all the test elements
      m_assembler.evaluateLocalWeakForms(TEST_TRIAL, m_testIndices, trialIndex,
                                         ALL_DOFS, localResult);

      // Global assembly
      {
//        MutexType::scoped_lock lock(m_mutex);
        // Loop over test indices
        for (int row = 0; row < testElementCount; ++row) {
          const int testIndex = m_testIndices[row];
          const int testDofCount = m_testGlobalDofs[testIndex].size();
          // Add the integrals to appropriate entries in the operator's matrix
          for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
            int trialGlobalDof = m_trialGlobalDofs[trialIndex][trialDof];
            if (trialGlobalDof < 0)
              continue;
            for (int testDof = 0; testDof < testDofCount; ++testDof) {
              int testGlobalDof = m_testGlobalDofs[testIndex][testDof];
              if (testGlobalDof < 0)
                continue;
              assert(std::abs(m_testLocalDofWeights[testIndex][testDof]) > 0.);
              assert(std::abs(m_trialLocalDofWeights[trialIndex][trialDof]) >
                     0.);
              MutexType::scoped_lock lock(m_mutex(testGlobalDof, trialGlobalDof));
              m_result(testGlobalDof, trialGlobalDof) +=
                  conj(m_testLocalDofWeights[testIndex][testDof]) *
                  m_trialLocalDofWeights[trialIndex][trialDof] *
                  localResult[row](testDof, trialDof);
            }
          }
        }
      }
    }
  }

private:
  const std::vector<int> &m_testIndices;
  const std::vector<std::vector<GlobalDofIndex>> &m_testGlobalDofs;
  const std::vector<std::vector<GlobalDofIndex>> &m_trialGlobalDofs;
  const std::vector<std::vector<BasisFunctionType>> &m_testLocalDofWeights;
  const std::vector<std::vector<BasisFunctionType>> &m_trialLocalDofWeights;
  // mutable OK because Assembler is thread-safe. (Alternative to "mutable"
  // here:
  // make assembler's internal integrator map mutable)
  typename Fiber::LocalAssemblerForIntegralOperators<ResultType> &m_assembler;
  // mutable OK because write access to this matrix is protected by a mutex
  Matrix<ResultType> &m_result;

  // mutex must be mutable because we need to lock and unlock it
  Matrix<MutexType> &m_mutex;
};

template <typename BasisFunctionType, typename ResultType>
class DensePotentialOperatorAssemblerLoopBody {
public:
  typedef tbb::spin_mutex MutexType;
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  DensePotentialOperatorAssemblerLoopBody(
      const std::vector<int> &pointIndices,
      const std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
      const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
      Fiber::LocalAssemblerForPotentialOperators<ResultType> &assembler,
      Matrix<ResultType> &result, MutexType &mutex)
      : m_pointIndices(pointIndices), m_trialGlobalDofs(trialGlobalDofs),
        m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler),
        m_result(result), m_mutex(mutex) {}

  void operator()(const tbb::blocked_range<int> &r) const {
    // In the current implementation we don't try to vary integration order
    // with point-element distance. So we'll supply -1, i.e. "unknown distance".
    const CoordinateType nominalDistance = -1.;

    const int pointCount = m_pointIndices.size();
    const int componentCount = m_assembler.resultDimension();

    std::vector<Matrix<ResultType>> localResult;
    for (int trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
      // Handle this trial element only if it contributes to any global DOFs.
      bool skipTrialElement = true;
      const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
      for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
        int trialGlobalDof = m_trialGlobalDofs[trialIndex][trialDof];
        if (trialGlobalDof >= 0) {
          skipTrialElement = false;
          break;
        }
      }
      if (skipTrialElement)
        continue;

      // Evaluate integrals over pairs of the current trial element and
      // all the points and components
      m_assembler.evaluateLocalContributions(
          m_pointIndices, trialIndex, ALL_DOFS, localResult, nominalDistance);

      // Global assembly
      {
        MutexType::scoped_lock lock(m_mutex);
        // Add the integrals to appropriate entries in the operator's matrix
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
          int trialGlobalDof = m_trialGlobalDofs[trialIndex][trialDof];
          if (trialGlobalDof < 0)
            continue;
          assert(std::abs(m_trialLocalDofWeights[trialIndex][trialDof]) > 0.);
          // Loop over point indices
          for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
            for (int component = 0; component < componentCount; ++component) {
              m_result(pointIndex * componentCount + component,
                       trialGlobalDof) +=
                  m_trialLocalDofWeights[trialIndex][trialDof] *
                  localResult[pointIndex](component, trialDof);
            }
          }
        }
      }
    }
  }

private:
  const std::vector<int> &m_pointIndices;
  const std::vector<std::vector<GlobalDofIndex>> &m_trialGlobalDofs;
  const std::vector<std::vector<BasisFunctionType>> &m_trialLocalDofWeights;
  // mutable OK because Assembler is thread-safe. (Alternative to "mutable"
  // here:
  // make assembler's internal integrator map mutable)
  typename Fiber::LocalAssemblerForPotentialOperators<ResultType> &m_assembler;
  // mutable OK because write access to this matrix is protected by a mutex
  Matrix<ResultType> &m_result;

  // mutex must be mutable because we need to lock and unlock it
  MutexType &m_mutex;
};

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid(). */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType> &space,
    std::vector<std::vector<GlobalDofIndex>> &globalDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) {
  // Get the grid's view so that we can iterate over elements
  const GridView &view = space.gridView();
  const int elementCount = view.entityCount(0);

  // Global DOF indices corresponding to local DOFs on elements
  globalDofs.clear();
  globalDofs.resize(elementCount);
  // Weights of the local DOFs on elements
  localDofWeights.clear();
  localDofWeights.resize(elementCount);

  // Gather global DOF lists
  const Mapper &mapper = view.elementMapper();
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    space.getGlobalDofs(element, globalDofs[elementIndex],
                        localDofWeights[elementIndex]);
    it->next();
  }
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
DenseGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &assembler,
    const Context<BasisFunctionType, ResultType> &context) {

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> testGlobalDofs, trialGlobalDofs;
  std::vector<std::vector<BasisFunctionType>> testLocalDofWeights,
      trialLocalDofWeights;
  gatherGlobalDofs(testSpace, testGlobalDofs, testLocalDofWeights);
  if (&testSpace == &trialSpace) {
    trialGlobalDofs = testGlobalDofs;
    trialLocalDofWeights = testLocalDofWeights;
  } else
    gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);
  const int testElementCount = testGlobalDofs.size();
  const int trialElementCount = trialGlobalDofs.size();

  // Enumerate the test elements that contribute to at least one global DOF
  std::vector<int> testIndices;
  testIndices.reserve(testElementCount);
  for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
    const int testDofCount = testGlobalDofs[testIndex].size();
    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      int testGlobalDof = testGlobalDofs[testIndex][testDof];
      if (testGlobalDof >= 0) {
        testIndices.push_back(testIndex);
        break;
      }
    }
  }

  // Create a discrete operator represented by a matrix that has to be calculated
  std::unique_ptr<DiscreteDenseBoundaryOperator<ResultType>>
  discreteDenseBoundaryOperator(new DiscreteDenseBoundaryOperator<ResultType>());

  // Create the operator's matrix
  Matrix<ResultType>& result = discreteDenseBoundaryOperator->matrix();
  result.resize(testSpace.globalDofCount(), trialSpace.globalDofCount());
  result.setZero();
  std::cout << "testGlobalDofCount = " << testSpace.globalDofCount()
      << ", trialGlobalDofCount = " << trialSpace.globalDofCount() << std::endl;

  typedef DenseWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
//  typename Body::MutexType mutex;

  // Create a mutex matrix for parallel assembly
  Matrix<typename Body::MutexType> mutex(testSpace.globalDofCount(),
                                         trialSpace.globalDofCount());
  {
    Fiber::SerialBlasRegion region;
    tbb::parallel_for(tbb::blocked_range<int>(0, trialElementCount),
                      Body(testIndices, testGlobalDofs, trialGlobalDofs,
                           testLocalDofWeights, trialLocalDofWeights, assembler,
                           result, mutex));
  }

  //// Old serial code (TODO: decide whether to keep it behind e.g. #ifndef
  /// PARALLEL)
  //    std::vector<arma::Mat<ValueType> > localResult;
  //    // Loop over trial elements
  //    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
  //    {
  //        // Evaluate integrals over pairs of the current trial element and
  //        // all the test elements
  //        assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices,
  //        trialIndex,
  //                                         ALL_DOFS, localResult);

  //        // Loop over test indices
  //        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
  //            // Add the integrals to appropriate entries in the operator's
  //            matrix
  //            for (int trialDof = 0; trialDof <
  //            trialGlobalDofs[trialIndex].size(); ++trialDof)
  //                for (int testDof = 0; testDof <
  //                testGlobalDofs[testIndex].size(); ++testDof)
  //                result(testGlobalDofs[testIndex][testDof],
  //                       trialGlobalDofs[trialIndex][trialDof]) +=
  //                        localResult[testIndex](testDof, trialDof);
  //    }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for classical dense assembly = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  std::ofstream file("cpu_dense_assembly_timer.dat", std::ios::out | std::ios::app);
  file << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

  // Return the discrete operator represented by the matrix that has just been
  // calculated
  return discreteDenseBoundaryOperator;
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
DenseGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
    const Matrix<CoordinateType> &points,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForPotentialOperators &assembler,
    const EvaluationOptions &options) {
  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> trialGlobalDofs;
  std::vector<std::vector<BasisFunctionType>> trialLocalDofWeights;
  gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);

  const int trialElementCount = trialGlobalDofs.size();
  const int pointCount = points.cols();
  const int componentCount = assembler.resultDimension();

  // Make a vector of all element indices
  std::vector<int> pointIndices(pointCount);
  for (int i = 0; i < pointCount; ++i)
    pointIndices[i] = i;

  // Create the operator's matrix
  Matrix<ResultType> result(pointCount * componentCount,
                            trialSpace.globalDofCount());
  result.setZero();

  typedef DensePotentialOperatorAssemblerLoopBody<BasisFunctionType, ResultType>
      Body;
  typename Body::MutexType mutex;

  {
    Fiber::SerialBlasRegion region;
    tbb::parallel_for(tbb::blocked_range<int>(0, trialElementCount),
                      Body(pointIndices, trialGlobalDofs, trialLocalDofWeights,
                           assembler, result, mutex));
  }
  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DenseGlobalAssembler);

} // namespace Bempp
