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

#include "elementary_integral_operator_base.hpp"

#include "context.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "numerical_quadrature_strategy.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>::
    ElementaryIntegralOperatorBase(
        const shared_ptr<const Space<BasisFunctionType>> &domain,
        const shared_ptr<const Space<BasisFunctionType>> &range,
        const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
        const std::string &label, int symmetry)
    : Base(domain, range, dualToRange, label, symmetry) {}

template <typename BasisFunctionType, typename ResultType>
ElementaryIntegralOperatorBase<BasisFunctionType,
                               ResultType>::~ElementaryIntegralOperatorBase() {}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>::
    assembleWeakFormInternal(
        LocalAssembler &assembler,
        const Context<BasisFunctionType, ResultType> &context) const {
  return assembleWeakFormInternalImpl2(assembler, context);
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<typename ElementaryIntegralOperatorBase<
    BasisFunctionType, ResultType>::LocalAssembler>
ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>::makeAssembler(
    const QuadratureStrategy &quadStrategy,
    const shared_ptr<const GeometryFactory> &testGeometryFactory,
    const shared_ptr<const GeometryFactory> &trialGeometryFactory,
    const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>> &
        testRawGeometry,
    const shared_ptr<const Fiber::RawGridGeometry<CoordinateType>> &
        trialRawGeometry,
    const shared_ptr<const std::vector<
        const Fiber::Shapeset<BasisFunctionType> *>> &testShapesets,
    const shared_ptr<const std::vector<
        const Fiber::Shapeset<BasisFunctionType> *>> &trialShapesets,
    const shared_ptr<const Fiber::OpenClHandler> &openClHandler,
    const ParallelizationOptions &parallelizationOptions,
    VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const {
  return makeAssemblerImpl(
      quadStrategy, testGeometryFactory, trialGeometryFactory, testRawGeometry,
      trialRawGeometry, testShapesets, trialShapesets, openClHandler,
      parallelizationOptions, verbosityLevel, cacheSingularIntegrals);
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<typename ElementaryIntegralOperatorBase<
    BasisFunctionType, ResultType>::LocalAssembler>
ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>::makeAssembler(
    const QuadratureStrategy &quadStrategy,
    const AssemblyOptions &options) const {
  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType> *>
      ShapesetPtrVector;

  const bool verbose = (options.verbosityLevel() >= VerbosityLevel::DEFAULT);

  shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
  shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
  shared_ptr<Fiber::OpenClHandler> openClHandler;
  shared_ptr<ShapesetPtrVector> testShapesets, trialShapesets;
  bool cacheSingularIntegrals;

  if (verbose)
    std::cout << "Collecting data for assembler construction..." << std::endl;
  this->collectDataForAssemblerConstruction(
      options, testRawGeometry, trialRawGeometry, testGeometryFactory,
      trialGeometryFactory, testShapesets, trialShapesets, openClHandler,
      cacheSingularIntegrals);
  if (verbose)
    std::cout << "Data collection finished." << std::endl;

  return makeAssemblerImpl(quadStrategy, testGeometryFactory,
                           trialGeometryFactory, testRawGeometry,
                           trialRawGeometry, testShapesets, trialShapesets,
                           openClHandler, options.parallelizationOptions(),
                           options.verbosityLevel(), cacheSingularIntegrals);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    ElementaryIntegralOperatorBase);

} // namespace Bempp
