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

#include "hmat_global_assembler.hpp"

#include "assembly_options.hpp"
#include "context.hpp"
#include "evaluation_options.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "weak_form_hmat_assembly_helper.hpp"
#include "discrete_hmat_boundary_operator.hpp"
#include "hmat_interface.hpp"

#include "../common/auto_timer.hpp"
#include "../common/chunk_statistics.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/shared_ptr.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../space/space.hpp"
#include "../common/bounding_box.hpp"

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/geometry_interface.hpp"
#include "../hmat/geometry_data_type.hpp"
#include "../hmat/geometry.hpp"
#include "../hmat/hmatrix.hpp"
#include "../hmat/data_accessor.hpp"
#include "../hmat/hmatrix_dense_compressor.hpp"
#include "../hmat/hmatrix_aca_compressor.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_queue.h>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
HMatGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    const std::vector<LocalAssemblerForIntegralOperators *> &localAssemblers,
    const std::vector<LocalAssemblerForIntegralOperators *>
        &localAssemblersForAdmissibleBlocks,
    const std::vector<const DiscreteBndOp *> &sparseTermsToAdd,
    const std::vector<ResultType> &denseTermMultipliers,
    const std::vector<ResultType> &sparseTermMultipliers,
    const Context<BasisFunctionType, ResultType> &context, int symmetry) {

  const AssemblyOptions &options = context.assemblyOptions();
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
  auto matVecParallelLevels = parameterList.template get<int>("options.hmat.matVecParallelLevels");
  if (coarseningAccuracy == 0)
    coarseningAccuracy = eps;

  shared_ptr<hmat::DefaultHMatrixType<ResultType>> hMatrix;

  if (compressionAlgorithm == "aca") {

    hmat::HMatrixAcaCompressor<ResultType, 2> compressor(helper, eps, maxRank);
    hMatrix.reset(new hmat::DefaultHMatrixType<ResultType>(
        blockClusterTree, compressor, matVecParallelLevels, coarsening, coarseningAccuracy));
  } else if (compressionAlgorithm == "dense") {
    hmat::HMatrixDenseCompressor<ResultType, 2> compressor(helper);
    hMatrix.reset(
        new hmat::DefaultHMatrixType<ResultType>(blockClusterTree, compressor, matVecParallelLevels));
  } else
    throw std::runtime_error("HMatGlobalAssember::assembleDetachedWeakForm: "
                             "Unknown compression algorithm");
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      static_cast<DiscreteBoundaryOperator<ResultType> *>(
          new DiscreteHMatBoundaryOperator<ResultType>(hMatrix)));
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
HMatGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &localAssembler,
    LocalAssemblerForIntegralOperators &localAssemblerForAdmissibleBlocks,
    const Context<BasisFunctionType, ResultType> &context, int symmetry) {
  typedef LocalAssemblerForIntegralOperators Assembler;
  std::vector<Assembler *> localAssemblers(1, &localAssembler);
  std::vector<Assembler *> localAssemblersForAdmissibleBlocks(
      1, &localAssemblerForAdmissibleBlocks);
  std::vector<const DiscreteBndOp *> sparseTermsToAdd;
  std::vector<ResultType> denseTermsMultipliers(1, 1.0);
  std::vector<ResultType> sparseTermsMultipliers;

  return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                                  localAssemblersForAdmissibleBlocks,
                                  sparseTermsToAdd, denseTermsMultipliers,
                                  sparseTermsMultipliers, context, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
HMatGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
        const Matrix<CoordinateType> &points,
                            const Space<BasisFunctionType> &trialSpace,
                            LocalAssemblerForPotentialOperators &localAssembler,
                            const ParameterList &parameterList)
{

    const size_t pointCount = points.cols();
    const int componentCount = localAssembler.resultDimension();
    const size_t testDofCount = pointCount * componentCount;
    const size_t trialDofCount = trialSpace.globalDofCount();



    return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>();
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(HMatGlobalAssembler);

} // namespace Bempp
