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


#include "elementary_local_operator.hpp"

#include "assembly_options.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "context.hpp"

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/local_assembler_for_local_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/eigen_support.hpp"
#include <boost/type_traits/is_complex.hpp>

#include <tbb/tick_count.h>
#include <tbb/tbb.h>

#include <stdexcept>
#include <vector>

namespace Bempp {

namespace {

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid(). */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    std::vector<std::vector<GlobalDofIndex>> &testGlobalDofs,
    std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
    std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
    std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights) {
  // We use the fact that test and trial space are required to be defined
  // on the same grid

  // Get the grid's leaf view so that we can iterate over elements
  const GridView &view = testSpace.gridView();
  const int elementCount = view.entityCount(0);

  // Global DOF indices corresponding to local DOFs on elements
  testGlobalDofs.clear();
  testGlobalDofs.resize(elementCount);
  trialGlobalDofs.clear();
  trialGlobalDofs.resize(elementCount);
  // Weights of the local DOFs on elements
  testLocalDofWeights.clear();
  testLocalDofWeights.resize(elementCount);
  trialLocalDofWeights.clear();
  trialLocalDofWeights.resize(elementCount);

  // Gather global DOF lists
  const Mapper &mapper = view.elementMapper();
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    testSpace.getGlobalDofs(element, testGlobalDofs[elementIndex],
                            testLocalDofWeights[elementIndex]);
    trialSpace.getGlobalDofs(element, trialGlobalDofs[elementIndex],
                             trialLocalDofWeights[elementIndex]);
    it->next();
  }
}

template <typename ValueType>
void generateTriplets(
    RealSparseMatrix &mat,
    const std::vector<std::vector<GlobalDofIndex>> &testGDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGDofs,
    const std::vector<Matrix<ValueType>> &localResult, int elementCount,
    std::vector<Eigen::Triplet<double>> &result);

template <>
void generateTriplets<double>(
    RealSparseMatrix &mat,
    const std::vector<std::vector<GlobalDofIndex>> &testGDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGDofs,
    const std::vector<Matrix<double>> &localResult, int elementCount,
    std::vector<Eigen::Triplet<double>> &result) {

  result.clear();
  size_t numberOfTriplets = 0;

  for (size_t e = 0; e < elementCount; ++e) {
    numberOfTriplets += localResult[e].rows() * localResult[e].cols();
  }

  result.reserve(numberOfTriplets);

  for (size_t e = 0; e < elementCount; ++e) {
    assert(testGDofs[e].size() == localResult[e].rows());
    assert(trialGDofs[e].size() == localResult[e].cols());
    for (int j = 0; j < localResult[e].cols(); ++j) {
      if (trialGDofs[e][j] < 0)
        continue;
      for (int i = 0; i < localResult[e].rows(); ++i) {
        if (testGDofs[e][i] < 0)
          continue;
        result.push_back(Eigen::Triplet<double>(
            testGDofs[e][i], trialGDofs[e][j], localResult[e](i, j)));
      }
    }
  }
}

template <>
void generateTriplets<float>(
    RealSparseMatrix &mat,
    const std::vector<std::vector<GlobalDofIndex>> &testGDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGDofs,
    const std::vector<Matrix<float>> &localResult, int elementCount,
    std::vector<Eigen::Triplet<double>> &result) {

  std::vector<Matrix<double>> localResult_tmp;
  localResult_tmp.reserve(localResult.size());
  for (int e = 0; e < elementCount; ++e)
    localResult_tmp.push_back(localResult[e].cast<double>());

  generateTriplets<double>(mat, testGDofs, trialGDofs, localResult_tmp,
                           elementCount, result);
}

template <>
void generateTriplets<std::complex<float>>(
    RealSparseMatrix &mat,
    const std::vector<std::vector<GlobalDofIndex>> &testGDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGDofs,
    const std::vector<Matrix<std::complex<float>>> &localResult,
    int elementCount, std::vector<Eigen::Triplet<double>> &result) {

  std::vector<Matrix<double>> localResult_tmp;
  localResult_tmp.reserve(localResult.size());
  for (int e = 0; e < elementCount; ++e)
    localResult_tmp.push_back(localResult[e].real().cast<double>());

  generateTriplets<double>(mat, testGDofs, trialGDofs, localResult_tmp,
                           elementCount, result);
}

template <>
void generateTriplets<std::complex<double>>(
    RealSparseMatrix &mat,
    const std::vector<std::vector<GlobalDofIndex>> &testGDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGDofs,
    const std::vector<Matrix<std::complex<double>>> &localResult,
    int elementCount, std::vector<Eigen::Triplet<double>> &result) {

  std::vector<Matrix<double>> localResult_tmp;
  localResult_tmp.reserve(localResult.size());
  for (int e = 0; e < elementCount; ++e)
    localResult_tmp.push_back(localResult[e].real());

  generateTriplets<double>(mat, testGDofs, trialGDofs, localResult_tmp,
                           elementCount, result);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// ElementaryLocalOperator

template <typename BasisFunctionType, typename ResultType>
ElementaryLocalOperator<BasisFunctionType, ResultType>::ElementaryLocalOperator(
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
    : Base(domain, range, dualToRange, label, symmetry) {}

template <typename BasisFunctionType, typename ResultType>
bool ElementaryLocalOperator<BasisFunctionType, ResultType>::isLocal() const {
  return true;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryLocalOperator<BasisFunctionType, ResultType>::assembleWeakFormImpl(
    const Context<BasisFunctionType, ResultType> &context) const {

  std::unique_ptr<LocalAssembler> assembler =
      this->makeAssembler(*context.quadStrategy(), context.assemblyOptions());
  shared_ptr<DiscreteBoundaryOperator<ResultType>> result =
      assembleWeakFormInternalImpl2(*assembler, context);
  tbb::tick_count::now();

  return result;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryLocalOperator<BasisFunctionType, ResultType>::
    assembleWeakFormInternalImpl2(
        LocalAssembler &assembler,
        const Context<BasisFunctionType, ResultType> &context) const {
  return shared_ptr<DiscreteBoundaryOperator<ResultType>>(
      assembleWeakFormInSparseMode(assembler, context.assemblyOptions())
          .release());
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryLocalOperator<BasisFunctionType, ResultType>::
    assembleWeakFormInDenseMode(LocalAssembler &assembler,
                                const AssemblyOptions &options) const {
  const Space<BasisFunctionType> &testSpace = *this->dualToRange();
  const Space<BasisFunctionType> &trialSpace = *this->domain();

  // Fill local submatrices
  const GridView &view = testSpace.gridView();
  const size_t elementCount = view.entityCount(0);
  std::vector<int> elementIndices(elementCount);
  for (size_t i = 0; i < elementCount; ++i)
    elementIndices[i] = i;
  std::vector<Matrix<ResultType>> localResult;
  assembler.evaluateLocalWeakForms(elementIndices, localResult);

  // Create the operator's matrix
  Matrix<ResultType> result(testSpace.globalDofCount(),
                            trialSpace.globalDofCount());
  result.setZero();

  // Retrieve global DOFs corresponding to local DOFs on all elements
  std::vector<std::vector<GlobalDofIndex>> testGdofs(elementCount);
  std::vector<std::vector<GlobalDofIndex>> trialGdofs(elementCount);
  std::vector<std::vector<BasisFunctionType>> testLdofWeights(elementCount);
  std::vector<std::vector<BasisFunctionType>> trialLdofWeights(elementCount);
  gatherGlobalDofs(testSpace, trialSpace, testGdofs, trialGdofs,
                   testLdofWeights, trialLdofWeights);

  // Distribute local matrices into the global matrix
  for (size_t e = 0; e < elementCount; ++e)
    for (size_t trialIndex = 0; trialIndex < trialGdofs[e].size();
         ++trialIndex) {
      int trialGdof = trialGdofs[e][trialIndex];
      if (trialGdof < 0)
        continue;
      for (size_t testIndex = 0; testIndex < testGdofs[e].size(); ++testIndex) {
        int testGdof = testGdofs[e][testIndex];
        if (testGdof < 0)
          continue;
        result(testGdof, trialGdof) += conj(testLdofWeights[e][testIndex]) *
                                       trialLdofWeights[e][trialIndex] *
                                       localResult[e](testIndex, trialIndex);
      }
    }

  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryLocalOperator<BasisFunctionType, ResultType>::
    assembleWeakFormInSparseMode(LocalAssembler &assembler,
                                 const AssemblyOptions &options) const {
  if (boost::is_complex<BasisFunctionType>::value)
    throw std::runtime_error(
        "ElementaryLocalOperator::assembleWeakFormInSparseMode(): "
        "sparse-mode assembly of identity operators for "
        "complex-valued basis functions is not supported yet");

  const Space<BasisFunctionType> &testSpace = *this->dualToRange();
  const Space<BasisFunctionType> &trialSpace = *this->domain();

  // Fill local submatrices
  const GridView &view = testSpace.gridView();
  const size_t elementCount = view.entityCount(0);
  std::vector<int> elementIndices(elementCount);
  for (size_t i = 0; i < elementCount; ++i)
    elementIndices[i] = i;
  std::vector<Matrix<ResultType>> localResult(elementCount);
  tbb::parallel_for(tbb::blocked_range<size_t>(0, elementCount),
          [&assembler,&localResult](tbb::blocked_range<size_t>& r){
          size_t numElements = r.end() - r.begin();
          std::vector<int> myIndices(numElements);
          std::vector<Matrix<ResultType>> myLocalResult;
          myLocalResult.reserve(numElements);
          int count = 0;
          for (size_t i = r.begin(); i!= r.end(); ++i)
            myIndices[count++] = i;
          assembler.evaluateLocalWeakForms(myIndices, myLocalResult);
          for (size_t i = 0; i < numElements; ++i)
            localResult[myIndices[i]] = myLocalResult[i];
            });
          
  //assembler.evaluateLocalWeakForms(elementIndices, localResult);

  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> testGdofs(elementCount);
  std::vector<std::vector<GlobalDofIndex>> trialGdofs(elementCount);
  std::vector<std::vector<BasisFunctionType>> testLdofWeights(elementCount);
  std::vector<std::vector<BasisFunctionType>> trialLdofWeights(elementCount);
  gatherGlobalDofs(testSpace, trialSpace, testGdofs, trialGdofs,
                   testLdofWeights, trialLdofWeights);

  // Multiply matrix entries by DOF weights
  for (size_t e = 0; e < elementCount; ++e)
    for (size_t trialDof = 0; trialDof < trialGdofs[e].size(); ++trialDof)
      for (size_t testDof = 0; testDof < testGdofs[e].size(); ++testDof)
        localResult[e](testDof, trialDof) *=
            conj(testLdofWeights[e][testDof]) * trialLdofWeights[e][trialDof];

  const int testGlobalDofCount = testSpace.globalDofCount();
  const int trialGlobalDofCount = trialSpace.globalDofCount();

  shared_ptr<RealSparseMatrix> result = boost::make_shared<RealSparseMatrix>(
      testGlobalDofCount, trialGlobalDofCount);

  std::vector<Eigen::Triplet<double>> triplets;
  generateTriplets<ResultType>(*result, testGdofs, trialGdofs, localResult,
                               elementCount, triplets);
  result->setFromTriplets(triplets.begin(), triplets.end());

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteSparseBoundaryOperator<ResultType>(result, this->symmetry(),
                                                     NO_TRANSPOSE));
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<typename ElementaryLocalOperator<BasisFunctionType,
                                                 ResultType>::LocalAssembler>
ElementaryLocalOperator<BasisFunctionType, ResultType>::makeAssembler(
    const QuadratureStrategy &quadStrategy,
    const AssemblyOptions &options) const {
  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType> *>
      ShapesetPtrVector;

  shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
  shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
  shared_ptr<Fiber::OpenClHandler> openClHandler;
  shared_ptr<ShapesetPtrVector> testShapesets, trialShapesets;
  bool cacheSingularIntegrals;

  this->collectDataForAssemblerConstruction(
      options, testRawGeometry, trialRawGeometry, testGeometryFactory,
      trialGeometryFactory, testShapesets, trialShapesets, openClHandler,
      cacheSingularIntegrals);
  assert(testRawGeometry == trialRawGeometry);
  assert(testGeometryFactory == trialGeometryFactory);

  return quadStrategy.makeAssemblerForLocalOperators(
      testGeometryFactory, testRawGeometry, testShapesets, trialShapesets,
      make_shared_from_ref(testTransformations()),
      make_shared_from_ref(trialTransformations()),
      make_shared_from_ref(integral()), openClHandler);
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<typename ElementaryLocalOperator<BasisFunctionType,
                                                 ResultType>::LocalAssembler>
ElementaryLocalOperator<BasisFunctionType, ResultType>::makeAssembler(
        const ParameterList& parameterList) const
{
    Context<BasisFunctionType, ResultType> context(parameterList);
    return this->makeAssembler(*context.quadStrategy(),context.assemblyOptions());

}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLocalOperator);

} // namespace Bempp
