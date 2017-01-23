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

#include "cuda_hmat_global_assembler.hpp"

#include "cuda_hmatrix_aca_compressor.hpp"

#include "../assembly/potential_operator_hmat_assembly_helper.hpp"
#include "../assembly/assembly_options.hpp"
#include "../assembly/context.hpp"
#include "../assembly/evaluation_options.hpp"
#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../assembly/weak_form_hmat_assembly_helper.hpp"
#include "../assembly/discrete_hmat_boundary_operator.hpp"
#include "../assembly/hmat_interface.hpp"

#include "../common/auto_timer.hpp"
#include "../common/chunk_statistics.hpp"
#include "../common/to_string.hpp"
#include "../common/bounding_box.hpp"
#include "../common/not_implemented_error.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/shared_ptr.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp"

#include "../grid/grid.hpp"

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/geometry_interface.hpp"
#include "../hmat/geometry_data_type.hpp"
#include "../hmat/geometry.hpp"
#include "../hmat/data_accessor.hpp"
#include "../hmat/hmatrix_dense_compressor.hpp"
#include "../hmat/hmatrix_data.hpp"
#include "../hmat/shared_ptr.hpp"

#include "../space/space.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <chrono>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
void CudaHMatGlobalAssembler<BasisFunctionType, ResultType>::
assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
    const Shapeset &testShapeset, const Shapeset &trialShapeset,
    const std::vector<LocalAssembler*> &localAssemblers,
    const std::vector<LocalAssembler*> &localAssemblersForAdmissibleBlocks,
    const std::vector<const DiscreteBoundaryOp*> &sparseTermsToAdd,
    const std::vector<ResultType> &denseTermMultipliers,
    const std::vector<ResultType> &sparseTermMultipliers,
    const Context<BasisFunctionType, ResultType> &context, int symmetry,
    shared_ptr<hmat::DefaultHMatrixType<ResultType>> &hMatrix) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const CudaOptions &cudaOptions = context.cudaOptions();
  const auto parameterList = context.globalParameterList();

  auto testSpacePointer = Fiber::make_shared_from_const_ref(testSpace);
  auto trialSpacePointer = Fiber::make_shared_from_const_ref(trialSpace);

  shared_ptr<const Space<BasisFunctionType>> actualTestSpace;
  shared_ptr<const Space<BasisFunctionType>> actualTrialSpace;
  actualTestSpace = testSpacePointer;
  actualTrialSpace = trialSpacePointer;

  auto minBlockSize =
      parameterList.template get<int>("options.hmat.minBlockSize");
  auto maxBlockSize =
      parameterList.template get<int>("options.hmat.maxBlockSize");
  auto eta = parameterList.template get<double>("options.hmat.eta");

  auto blockClusterTree = generateBlockClusterTree(
      *actualTestSpace, *actualTrialSpace, parameterList);

  WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType> helper(
      *actualTestSpace, *actualTrialSpace, blockClusterTree, localAssemblers,
      sparseTermsToAdd, denseTermMultipliers, sparseTermMultipliers);

  auto compressionAlgorithm = parameterList.template get<std::string>(
      "options.hmat.compressionAlgorithm");

  auto maxRank = parameterList.template get<int>("options.hmat.maxRank");
  auto eps = parameterList.template get<double>("options.hmat.eps");
  auto coarsening = parameterList.template get<bool>("options.hmat.coarsening");
  auto coarseningAccuracy =
      parameterList.template get<double>("options.hmat.coarseningAccuracy");
  auto matVecParallelLevels =
      parameterList.template get<int>("options.hmat.matVecParallelLevels");
  if (coarseningAccuracy == 0)
    coarseningAccuracy = 0.1 * eps; // Default is finer tol for coarsening than eps

  if (compressionAlgorithm == "aca") {

    const int N = 2;

    typedef tbb::concurrent_unordered_map<
        shared_ptr<hmat::BlockClusterTreeNode<N>>,
        shared_ptr<hmat::HMatrixData<ResultType>>,
        hmat::shared_ptr_hash<hmat::BlockClusterTreeNode<N>>>
            ParallelDataContainer;
    ParallelDataContainer hMatrixData;
    hMatrixData.clear();

    // Get numerical quadrature points and weights
    const int trialQuadOrder = cudaOptions.quadOrder();
    const int testQuadOrder = cudaOptions.quadOrder();
    Matrix<CoordinateType> localTrialQuadPoints, localTestQuadPoints;
    std::vector<CoordinateType> trialQuadWeights, testQuadWeights;
    Fiber::fillSingleQuadraturePointsAndWeights(3, trialQuadOrder,
        localTrialQuadPoints, trialQuadWeights);
    Fiber::fillSingleQuadraturePointsAndWeights(3, testQuadOrder,
        localTestQuadPoints, testQuadWeights);

    if (context.cudaOptions().precision() == "single") {

      typedef typename ScalarTraits<typename ScalarTraits<BasisFunctionType>::RealType>::SingleType CudaBasisFunctionType;
      typedef typename ScalarTraits<typename ScalarTraits<KernelType>::RealType>::SingleType CudaKernelType;
      typedef typename ScalarTraits<typename ScalarTraits<ResultType>::RealType>::SingleType CudaResultType;

      hmat::CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType, N> compressor(
          *actualTestSpace, *actualTrialSpace, kernel, testShapeset, trialShapeset,
          localTestQuadPoints, localTrialQuadPoints,
          testQuadWeights, trialQuadWeights,
          blockClusterTree, localAssemblers, sparseTermsToAdd,
          denseTermMultipliers, sparseTermMultipliers,
          eps, 1e-15, maxRank, helper, cudaOptions);

      compressor.compressAllBlocks(hMatrixData);

    } else {

      typedef typename ScalarTraits<typename ScalarTraits<BasisFunctionType>::RealType>::DoubleType CudaBasisFunctionType;
      typedef typename ScalarTraits<typename ScalarTraits<KernelType>::RealType>::DoubleType CudaKernelType;
      typedef typename ScalarTraits<typename ScalarTraits<ResultType>::RealType>::DoubleType CudaResultType;

      hmat::CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType, N> compressor(
          *actualTestSpace, *actualTrialSpace, kernel, testShapeset, trialShapeset,
          localTestQuadPoints, localTrialQuadPoints,
          testQuadWeights, trialQuadWeights,
          blockClusterTree, localAssemblers, sparseTermsToAdd,
          denseTermMultipliers, sparseTermMultipliers,
          eps, 1e-15, maxRank, helper, cudaOptions);

      compressor.compressAllBlocks(hMatrixData);
    }

    hMatrix.reset(new hmat::DefaultHMatrixType<ResultType>(
        hMatrixData, blockClusterTree, coarsening, coarseningAccuracy));

  } else if (compressionAlgorithm == "dense") {

    throw NotImplementedError(
        "CudaHMatGlobalAssember::assembleDetachedWeakForm: "
        "Dense compression algorithm is not currently supported");

  } else {

    throw std::runtime_error("CudaHMatGlobalAssember::assembleDetachedWeakForm: "
                             "Unknown compression algorithm");
  }
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaHMatGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssembler &localAssembler,
    LocalAssembler &localAssemblerForAdmissibleBlocks,
    const Context<BasisFunctionType, ResultType> &context, int symmetry) {

  std::chrono::steady_clock::time_point beginTotal = std::chrono::steady_clock::now();

  typedef LocalAssembler Assembler;

  std::vector<Assembler*> localAssemblers(1, &localAssembler);
  std::vector<Assembler*> localAssemblersForAdmissibleBlocks(
      1, &localAssemblerForAdmissibleBlocks);
  std::vector<const DiscreteBoundaryOp*> sparseTermsToAdd;
  std::vector<ResultType> denseTermsMultipliers(1, 1.0);
  std::vector<ResultType> sparseTermsMultipliers;

  // Cast assembler to default assembler
  const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
      BasisFunctionType, KernelType, ResultType, GeometryFactory>
  &defaultAssembler =
      reinterpret_cast<const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
          BasisFunctionType, KernelType, ResultType, GeometryFactory> &>(localAssembler);

  // Get kernel from assembler
  shared_ptr<const Fiber::CollectionOfKernels<KernelType>> kernel =
      defaultAssembler.kernels();

  // Get shapesets from assembler
  shared_ptr<const std::vector<const Shapeset*>>
  testShapesets, trialShapesets;
  defaultAssembler.getShapesets(testShapesets, trialShapesets);
  const Shapeset &testShapeset = *(*testShapesets)[0];
  const Shapeset &trialShapeset = *(*trialShapesets)[0];
  std::cout << "NOTE: Shapesets have to be identical within one space" << std::endl;
  std::cout << "trialDofCount = " << trialShapeset.size() << ", " << std::flush;
  std::cout << "testDofCount = " << testShapeset.size() << std::endl;

  shared_ptr<hmat::DefaultHMatrixType<ResultType>> hMatrix;

  assembleDetachedWeakForm(
          testSpace, trialSpace, kernel, testShapeset, trialShapeset,
          localAssemblers, localAssemblersForAdmissibleBlocks,
          sparseTermsToAdd, denseTermsMultipliers,
          sparseTermsMultipliers, context, symmetry, hMatrix);

  std::chrono::steady_clock::time_point endTotal = std::chrono::steady_clock::now();
  std::cout << "Time for CUDA H-matrix assembly = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - beginTotal).count()
            << " ms" << std::endl;

  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      static_cast<DiscreteBoundaryOperator<ResultType> *>(
          new DiscreteHMatBoundaryOperator<ResultType>(hMatrix)));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaHMatGlobalAssembler);

} // namespace Bempp
