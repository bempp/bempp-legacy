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
#include "bempp/common/config_trilinos.hpp"

#include "hmat_global_assembler.hpp"

#include "assembly_options.hpp"
#include "context.hpp"
#include "evaluation_options.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "discrete_sparse_boundary_operator.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/auto_timer.hpp"
#include "../common/chunk_statistics.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/shared_ptr.hpp"
#include "../space/space.hpp"
#include "../common/bounding_box.hpp"

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/geometry_interface.hpp"
#include "../hmat/geometry_data_type.hpp"
#include "../hmat/geometry.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>


namespace Bempp {


namespace {

template <typename BasisFunctionType>
class SpaceHMatGeometryInterface : public hmat::GeometryInterface {

public:
  typedef typename Fiber::ScalarTraits<BasisFunctionType>::RealType
  CoordinateType;
  SpaceHMatGeometryInterface(const Space<BasisFunctionType> &space)
      : m_counter(0) {
    space.getGlobalDofBoundingBoxes(m_bemppBoundingBoxes);
  }
  shared_ptr<const hmat::GeometryDataType> next() const override {

    if (m_counter == m_bemppBoundingBoxes.size())
      return shared_ptr<hmat::GeometryDataType>();

    auto lbound = m_bemppBoundingBoxes[m_counter].lbound;
    auto ubound = m_bemppBoundingBoxes[m_counter].ubound;
    auto center = m_bemppBoundingBoxes[m_counter].reference;
    m_counter++;
    return shared_ptr<hmat::GeometryDataType>(new hmat::GeometryDataType(
        hmat::BoundingBox(lbound[0], ubound[0], lbound[1], ubound[1], lbound[2],
                          ubound[2]),
        std::array<double, 3>({{center.x, center.y, center.z}})));
  }

  void reset() override { m_counter = 0; }

private:
  std::size_t m_counter;
  std::vector<const BoundingBox<CoordinateType>> m_bemppBoundingBoxes;
};

template <typename BasisFunctionType>
shared_ptr<hmat::DefaultBlockClusterTreeType>
generateBlockClusterTree(const Space<BasisFunctionType> &testSpace,
                         const Space<BasisFunctionType> &trialSpace,
                         int minBlockSize, int maxBlockSize, double eta) {

  hmat::Geometry testGeometry;
  hmat::Geometry trialGeometry;

  hmat::fillGeometry(testGeometry,
                     SpaceHMatGeometryInterface<BasisFunctionType>(testSpace));
  hmat::fillGeometry(trialGeometry,
                     SpaceHMatGeometryInterface<BasisFunctionType>(trialSpace));

  auto testClusterTree = shared_ptr<hmat::DefaultClusterTreeType>(
      new hmat::DefaultClusterTreeType(testGeometry, minBlockSize));

  auto trialClusterTree = shared_ptr<hmat::DefaultClusterTreeType>(
      new hmat::DefaultClusterTreeType(trialGeometry, minBlockSize));

  return shared_ptr<hmat::DefaultBlockClusterTreeType>(
      new hmat::DefaultBlockClusterTreeType(testClusterTree, trialClusterTree,
                                            maxBlockSize,
                                            hmat::StandardAdmissibility(eta)));
}
} // end anonymous namespace
template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
HMatGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    const std::vector<LocalAssemblerForIntegralOperators *> &localAssemblers,
    const std::vector<LocalAssemblerForIntegralOperators *> &
        localAssemblersForAdmissibleBlocks,
    const std::vector<const DiscreteBndOp *> &sparseTermsToAdd,
    const std::vector<ResultType> &denseTermMultipliers,
    const std::vector<ResultType> &sparseTermMultipliers,
    const Context<BasisFunctionType, ResultType> &context, int symmetry) {

  bool symmetric = symmetry & SYMMETRIC;
  if (symmetry & HERMITIAN && !(symmetry & SYMMETRIC) &&
      verbosityAtLeastDefault)
    std::cout << "Warning: assembly of non-symmetric Hermitian H-matrices "
                 "is not supported yet. A general H-matrix will be assembled"
              << std::endl;

  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>();
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(HMatGlobalAssembler);

} // namespace Bempp
