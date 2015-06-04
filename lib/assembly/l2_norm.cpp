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

#include "l2_norm.hpp"

#include "evaluation_options.hpp"
#include "grid_function.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/eigen_support.hpp"
#include "../common/complex_aux.hpp"
#include "../common/shared_ptr.hpp"

#include "../fiber/collection_of_2d_arrays.hpp"
#include "../fiber/collection_of_4d_arrays.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_kernel_trial_integral.hpp"
#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/function.hpp"
#include "../fiber/geometrical_data.hpp"
#include "../grid/mapper.hpp"

namespace Bempp {

namespace {

template <typename ValueType_> class KernelFunctorFromUnaryFunction {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  explicit KernelFunctorFromUnaryFunction(
      const Fiber::Function<ValueType> &function)
      : m_function(function) {}

  int kernelCount() const { return 1; }

  int kernelRowCount(int kernelIndex) const {
    return m_function.codomainDimension();
  }

  int kernelColCount(int kernelIndex) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    m_function.addGeometricalDependencies(trialGeomDeps);
  }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(
      const Fiber::ConstGeometricalDataSlice<CoordinateType> &testGeomData,
      const Fiber::ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
      CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    Fiber::GeometricalData<CoordinateType> geomData =
        trialGeomData.asGeometricalData();
    Matrix<ValueType> resultMatrix;
    m_function.evaluate(geomData, resultMatrix);
    assert(resultMatrix.rows() == result[0].extent(0));
    assert(resultMatrix.cols() == result[0].extent(1));
    for (size_t j = 0; j < resultMatrix.cols(); ++j)
      for (size_t i = 0; i < resultMatrix.rows(); ++i)
        result[0](i, j) = resultMatrix(i, j);
  }

private:
  const Fiber::Function<ValueType> &m_function;
};

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class L2NormOfDifferenceIntegrandFunctor {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef KernelType_ KernelType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  void addGeometricalDependencies(size_t &trialGeomDeps) const {
    // do nothing
  }

  int resultDimension() const { return 1; }

  template <typename CollectionOf1dSlicesOfConstNdArrays>
  void evaluate(const Fiber::ConstGeometricalDataSlice<CoordinateType>
                    & /* trialGeomData */,
                const Fiber::CollectionOf2dSlicesOfConst4dArrays<KernelType>
                    &kernelValues,
                const CollectionOf1dSlicesOfConstNdArrays &trialValues,
                std::vector<ResultType> &result) const {
    // Assert that there is only a single kernel with a single column
    assert(kernelValues.size() == 1);
    assert(kernelValues[0].extent(1) == 1);
    const size_t componentCount = kernelValues[0].extent(0);

    // Assert that there is only one trial transformation
    // and that its dimensions are correct
    assert(trialValues.size() == 1);
#ifndef NDEBUG
    assert(trialValues[0].extent(0) == componentCount);
#endif

    // Assert that the result has only one element
    assert(result.size() == 1);

    // (t - k)* . (t - k)
    result[0] = 0.;
    for (size_t i = 0; i < componentCount; ++i) {
      result[0] += conj(trialValues[0](i)-kernelValues[0](i, 0)) *
                   (trialValues[0](i)-kernelValues[0](i, 0));
    }
  }
};

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<Fiber::EvaluatorForIntegralOperators<ResultType>>
makeEvaluator(const GridFunction<BasisFunctionType, ResultType> &gridFunction,
              const Fiber::Function<ResultType> &refFunction,
              const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                              GeometryFactory> &quadStrategy,
              const EvaluationOptions &options) {
  // Collect the standard set of data necessary for construction of
  // evaluators and assemblers
  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType CoordinateType;
  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType> *>
      ShapesetPtrVector;
  typedef std::vector<std::vector<ResultType>> CoefficientsVector;
  typedef LocalAssemblerConstructionHelper Helper;

  shared_ptr<RawGridGeometry> rawGeometry;
  shared_ptr<GeometryFactory> geometryFactory;
  shared_ptr<Fiber::OpenClHandler> openClHandler;
  shared_ptr<ShapesetPtrVector> shapesets;

  const Space<BasisFunctionType> &space = *gridFunction.space();
  Helper::collectGridData(space, rawGeometry, geometryFactory);
  Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                            rawGeometry, openClHandler);
  Helper::collectShapesets(space, shapesets);

  // In addition, get coefficients of argument's expansion in each element
  const GridView &view = space.gridView();
  const int elementCount = view.entityCount(0);

  shared_ptr<CoefficientsVector> localCoefficients =
      boost::make_shared<CoefficientsVector>(elementCount);

  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  const Mapper &mapper = view.elementMapper();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    gridFunction.getLocalCoefficients(element,
                                      (*localCoefficients)[elementIndex]);
    it->next();
  }

  // Construct kernels collection
  typedef KernelFunctorFromUnaryFunction<ResultType> KernelFunctor;
  typedef Fiber::DefaultCollectionOfKernels<KernelFunctor> Kernels;
  shared_ptr<Fiber::CollectionOfKernels<ResultType>> kernels =
      boost::make_shared<Kernels>(KernelFunctor(refFunction));

  // Construct trial function transformations collection
  const Fiber::CollectionOfShapesetTransformations<CoordinateType>
      &trialTransformations = space.basisFunctionValue();

  // Construct integral
  typedef L2NormOfDifferenceIntegrandFunctor<BasisFunctionType, ResultType,
                                             ResultType> IntegralFunctor;
  shared_ptr<Fiber::KernelTrialIntegral<BasisFunctionType, ResultType,
                                        ResultType>> integral =
      boost::make_shared<Fiber::DefaultKernelTrialIntegral<IntegralFunctor>>(
          IntegralFunctor());

  // Now create the evaluator
  return quadStrategy.makeEvaluatorForIntegralOperators(
      geometryFactory, rawGeometry, shapesets, kernels,
      make_shared_from_ref(trialTransformations), // lives as long as space does
      integral, localCoefficients, openClHandler,
      options.parallelizationOptions());
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType L2NormOfDifference(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    const Fiber::Function<ResultType> &refFunction,
    const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                    GeometryFactory> &quadStrategy,
    const EvaluationOptions &options) {
  // First, construct the evaluator.
  typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;
  std::unique_ptr<Evaluator> evaluator =
      makeEvaluator(gridFunction, refFunction, quadStrategy, options);

  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType MagnitudeType;
  typedef MagnitudeType CoordinateType;
  Matrix<CoordinateType> evaluationPoints(1, 1);
  evaluationPoints(0, 0) = 0.;
  Matrix<ResultType> resultMatrix;
  evaluator->evaluate(Evaluator::FAR_FIELD, evaluationPoints, resultMatrix);

  ResultType result = resultMatrix(0, 0);
  if (fabs(imagPart(result)) >
      1000. * std::numeric_limits<MagnitudeType>::epsilon())
    std::cout << "Warning: squared L2 norm has non-negligible imaginary part: "
              << imagPart(result) << std::endl;
  return sqrt(realPart(result));
}

template <typename BasisFunctionType, typename ResultType>
void estimateL2Error(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    const Fiber::Function<ResultType> &refFunction,
    const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                    GeometryFactory> &quadStrategy,
    const EvaluationOptions &options,
    typename ScalarTraits<BasisFunctionType>::RealType &absError,
    typename ScalarTraits<BasisFunctionType>::RealType &relError) {
  typedef
      typename Fiber::ScalarTraits<BasisFunctionType>::RealType MagnitudeType;
  absError =
      L2NormOfDifference(gridFunction, refFunction, quadStrategy, options);
  MagnitudeType refNorm =
      L2NormOfDifference(0. * gridFunction, refFunction, quadStrategy, options);
  relError = absError / refNorm;
}

template <typename BasisFunctionType, typename ResultType>
void estimateL2Error(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    const Fiber::Function<ResultType> &refFunction,
    const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                    GeometryFactory> &quadStrategy,
    typename ScalarTraits<BasisFunctionType>::RealType &absError,
    typename ScalarTraits<BasisFunctionType>::RealType &relError) {
  estimateL2Error(gridFunction, refFunction, quadStrategy, EvaluationOptions(),
                  absError, relError);
}

#define INSTANTIATE_FUNCTION(BASIS, RESULT)                                    \
  template ScalarTraits<BASIS>::RealType L2NormOfDifference(                   \
      const GridFunction<BASIS, RESULT> &gridFunction,                         \
      const Fiber::Function<RESULT> &refFunction,                              \
      const Fiber::QuadratureStrategy<BASIS, RESULT, GeometryFactory>          \
          &quadStrategy,                                                       \
      const EvaluationOptions &options);                                       \
  template void estimateL2Error(                                               \
      const GridFunction<BASIS, RESULT> &gridFunction,                         \
      const Fiber::Function<RESULT> &refFunction,                              \
      const Fiber::QuadratureStrategy<BASIS, RESULT, GeometryFactory>          \
          &quadStrategy,                                                       \
      const EvaluationOptions &options,                                        \
      ScalarTraits<BASIS>::RealType &absError,                                 \
      ScalarTraits<BASIS>::RealType &relError);                                \
  template void estimateL2Error(                                               \
      const GridFunction<BASIS, RESULT> &gridFunction,                         \
      const Fiber::Function<RESULT> &refFunction,                              \
      const Fiber::QuadratureStrategy<BASIS, RESULT, GeometryFactory>          \
          &quadStrategy,                                                       \
      ScalarTraits<BASIS>::RealType &absError,                                 \
      ScalarTraits<BASIS>::RealType &relError)

FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_FUNCTION);

} // namespace Bempp
