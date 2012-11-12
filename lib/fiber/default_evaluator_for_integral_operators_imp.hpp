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

#include "default_evaluator_for_integral_operators.hpp" // keep IDEs happy

#include "../common/common.hpp"

#include "basis.hpp"
#include "basis_data.hpp"
#include "collection_of_basis_transformations.hpp"
#include "geometrical_data.hpp"
#include "collection_of_kernels.hpp"
#include "collection_of_2d_arrays.hpp"
#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "kernel_trial_integral.hpp"
#include "numerical_quadrature.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "serial_blas_region.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace Fiber
{

namespace
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class EvaluationLoopBody
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    EvaluationLoopBody(
            size_t chunkSize,
            const arma::Mat<CoordinateType>& points,
            const GeometricalData<CoordinateType>& trialGeomData,
            const CollectionOf2dArrays<ResultType>& trialTransfValues,
            const std::vector<CoordinateType>& weights,
            const CollectionOfKernels<KernelType>& kernels,
            const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>& integral,
            arma::Mat<ResultType>& result) :
        m_chunkSize(chunkSize),
        m_points(points), m_trialGeomData(trialGeomData),
        m_trialTransfValues(trialTransfValues), m_weights(weights),
        m_kernels(kernels), m_integral(integral), m_result(result),
        m_pointCount(result.n_cols), m_outputComponentCount(result.n_rows)
    {
    }

    void operator() (const tbb::blocked_range<size_t>& r) const {
        CollectionOf4dArrays<KernelType> kernelValues;
        GeometricalData<CoordinateType> evalPointGeomData;
        for (size_t i = r.begin(); i < r.end(); ++i)
        {
            size_t start = m_chunkSize * i;
            size_t end = std::min(start + m_chunkSize, m_pointCount);
            evalPointGeomData.globals = m_points.cols(start, end - 1 /* inclusive */);
            m_kernels.evaluateOnGrid(evalPointGeomData, m_trialGeomData, kernelValues);
            // View into the current chunk of the "result" array
            _2dArray<ResultType> resultChunk(m_outputComponentCount, end - start,
                                             m_result.colptr(start));
            m_integral.evaluate(m_trialGeomData,
                                kernelValues,
                                m_trialTransfValues,
                                m_weights,
                                resultChunk);
        }
    }

private:
    size_t m_chunkSize;
    const arma::Mat<CoordinateType>& m_points;
    const GeometricalData<CoordinateType>& m_trialGeomData;
    const CollectionOf2dArrays<ResultType>& m_trialTransfValues;
    const std::vector<CoordinateType>& m_weights;
    const CollectionOfKernels<KernelType>& m_kernels;
    const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>& m_integral;
    arma::Mat<ResultType>& m_result;
    size_t m_pointCount;
    size_t m_outputComponentCount;
};

} // namespace

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::DefaultEvaluatorForIntegralOperators(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const CollectionOfKernels<KernelType> >& kernels,
        const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType> >& integral,
        const shared_ptr<const std::vector<std::vector<ResultType> > >& argumentLocalCoefficients,
        const shared_ptr<const OpenClHandler >& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        const QuadratureOptions& quadratureOptions) :
    m_geometryFactory(geometryFactory), m_rawGeometry(rawGeometry),
    m_trialBases(trialBases), m_kernels(kernels),
    m_trialTransformations(trialTransformations), m_integral(integral),
    m_argumentLocalCoefficients(argumentLocalCoefficients),
    m_openClHandler(openClHandler),
    m_parallelizationOptions(parallelizationOptions),
    m_quadratureOptions(quadratureOptions)
{
    const size_t elementCount = rawGeometry->elementCount();
    if (!rawGeometry->auxData().is_empty() &&
            rawGeometry->auxData().n_cols != elementCount)
        throw std::invalid_argument(
                "DefaultEvaluatorForIntegralOperators::"
                "DefaultEvaluatorForIntegralOperators(): "
                "number of columns of auxData must match that of "
                "elementCornerIndices");
    if (trialBases->size() != elementCount)
        throw std::invalid_argument(
                "DefaultEvaluatorForIntegralOperators::"
                "DefaultEvaluatorForIntegralOperators(): "
                "size of testBases must match the number of columns of "
                "elementCornerIndices");

    // Cache "trial data" such as the values of the argument at quadrature
    // points, vectors normal to surface at quadrature points etc.
    cacheTrialData();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::evaluate(
        Region region,
        const arma::Mat<CoordinateType>& points, arma::Mat<ResultType>& result) const
{
    const size_t pointCount = points.n_cols;
    const int outputComponentCount = m_integral->resultDimension();

    result.set_size(outputComponentCount, pointCount);
    result.fill(0.);

    const GeometricalData<CoordinateType>& trialGeomData =
            (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD) ?
                m_nearFieldTrialGeomData :
                m_farFieldTrialGeomData;
    const CollectionOf2dArrays<ResultType>& trialTransfValues =
            (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD) ?
                m_nearFieldTrialTransfValues :
                m_farFieldTrialTransfValues;
    const std::vector<CoordinateType>& weights =
            (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD) ?
                m_nearFieldWeights :
                m_farFieldWeights;

    // Do things in chunks of 96 points -- in order to avoid creating
    // too large arrays of kernel values
    const size_t chunkSize = 96;
    const size_t chunkCount = (pointCount + chunkSize - 1) / chunkSize;

    int maxThreadCount = 1;
    if (!m_parallelizationOptions.isOpenClEnabled()) {
        if (m_parallelizationOptions.maxThreadCount() ==
                ParallelizationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = m_parallelizationOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    typedef EvaluationLoopBody<
            BasisFunctionType, KernelType, ResultType> Body;
    {
        Fiber::SerialBlasRegion region;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, chunkCount),
                          Body(chunkSize,
                               points, trialGeomData, trialTransfValues, weights,
                               *m_kernels, *m_integral, result));
    }

//    // Old serial version
//    CollectionOf4dArrays<KernelType> kernelValues;
//    GeometricalData<CoordinateType> evalPointGeomData;
//    for (size_t start = 0; start < pointCount; start += chunkSize)
//    {
//        size_t end = std::min(start + chunkSize, pointCount);
//        evalPointGeomData.globals = points.cols(start, end - 1 /* inclusive */);
//        m_kernels->evaluateOnGrid(evalPointGeomData, trialGeomData, kernelValues);
//        // View into the current chunk of the "result" array
//        _2dArray<ResultType> resultChunk(outputComponentCount, end - start,
//                                         result.colptr(start));
//        m_integral->evaluate(trialGeomData,
//                             kernelValues,
//                             weightedTrialTransfValues,
//                             resultChunk);
//    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::cacheTrialData()
{
    size_t testGeomDeps = 0, trialGeomDeps = 0;
    m_kernels->addGeometricalDependencies(testGeomDeps, trialGeomDeps);
    if (testGeomDeps != 0 && testGeomDeps != GLOBALS)
        throw std::runtime_error(
                "DefaultEvaluatorForIntegralOperators::cacheTrialData(): "
                "potentials cannot contain kernels that depend on other test data "
                "than global coordinates");

    calcTrialData(EvaluatorForIntegralOperators<ResultType>::FAR_FIELD,
                  trialGeomDeps, m_farFieldTrialGeomData,
                  m_farFieldTrialTransfValues, m_farFieldWeights);
    // near field currently not used
    calcTrialData(EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD,
                  trialGeomDeps, m_nearFieldTrialGeomData,
                  m_nearFieldTrialTransfValues, m_nearFieldWeights);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::calcTrialData(
        Region region,
        int kernelTrialGeomDeps,
        GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf2dArrays<ResultType>& trialTransfValues,
        std::vector<CoordinateType>& weights) const
{
    const int elementCount = m_rawGeometry->elementCount();
    const int worldDim = m_rawGeometry->worldDimension();
    const int gridDim = m_rawGeometry->gridDimension();
    const int transformationCount = m_trialTransformations->transformationCount();

    // Find out which basis data need to be calculated
    size_t basisDeps = 0;
    // Find out which geometrical data need to be calculated, in addition
    // to those needed by the kernel
    size_t trialGeomDeps = kernelTrialGeomDeps;
    m_trialTransformations->addDependencies(basisDeps, trialGeomDeps);
    trialGeomDeps |= INTEGRATION_ELEMENTS;

    // Initialise the geometry factory
    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory->make());

    // Find all unique trial bases
    // Set of unique quadrature variants
    typedef std::set<const Basis<BasisFunctionType>*> BasisSet;
    BasisSet uniqueTrialBases(m_trialBases->begin(), m_trialBases->end());

    // Initialise temporary (per-element) data containers
    std::vector<GeometricalData<CoordinateType> > geomDataPerElement(elementCount);
    std::vector<CollectionOf2dArrays<ResultType> >
            trialTransfValuesPerElement(elementCount);
    std::vector<std::vector<CoordinateType> > weightsPerElement(elementCount);

    int quadPointCount = 0; // Quadrature point counter

    for (typename BasisSet::const_iterator it = uniqueTrialBases.begin();
         it != uniqueTrialBases.end(); ++it)
    {
        const Basis<BasisFunctionType>& activeBasis = *(*it);
        int order = quadOrder(activeBasis, region);

        // Find out the element type
        int elementCornerCount = 0;
        for (int e = 0; e < elementCount; ++e)
            if ((*m_trialBases)[e] == &activeBasis)
            {
                // elementCornerCount = m_rawGeometry->elementCornerCount(e);
                // This implementation prevents a segmentation fault on Macs
                // when compiled with llvm in 64-bit mode with -O2 or -O3
                const arma::Mat<int>& elementCornerIndices =
                    m_rawGeometry->elementCornerIndices();
                for (size_t i = 0; i < elementCornerIndices.n_rows; ++i)
                    if (elementCornerIndices(i, e) >= 0)
                        elementCornerCount = i + 1;
                    else
                        break;

                break;
            }

        // Get quadrature points and weights
        arma::Mat<CoordinateType> localQuadPoints;
        std::vector<CoordinateType> quadWeights;
        fillSingleQuadraturePointsAndWeights(
                    elementCornerCount, order, localQuadPoints, quadWeights);

        // Get basis data
        BasisData<BasisFunctionType> basisData;
        activeBasis.evaluate(basisDeps, localQuadPoints, ALL_DOFS, basisData);

        BasisData<ResultType> argumentData;
        if (basisDeps & VALUES)
            argumentData.values.set_size(basisData.values.extent(0),
                                         1, // just one function
                                         basisData.values.extent(2));
        if (basisDeps & DERIVATIVES)
            argumentData.derivatives.set_size(basisData.derivatives.extent(0),
                                              basisData.derivatives.extent(1),
                                              1, // just one function
                                              basisData.derivatives.extent(3));

        // Loop over elements and process those that use the active basis
        CollectionOf3dArrays<ResultType> trialValues;
        for (int e = 0; e < elementCount; ++e)
        {
            if ((*m_trialBases)[e] != &activeBasis)
                continue;

            // Local coefficients of the argument in the current element
            const std::vector<ResultType>& localCoefficients =
                    (*m_argumentLocalCoefficients)[e];

            // Calculate the argument function's values and/or derivatives
            // at quadrature points in the current element
            if (basisDeps & VALUES)
            {
                std::fill(argumentData.values.begin(),
                          argumentData.values.end(), 0.);
                for (size_t point = 0; point < basisData.values.extent(2); ++point)
                    for (size_t dim = 0; dim < basisData.values.extent(0); ++dim)
                        for (size_t fun = 0; fun < basisData.values.extent(1); ++fun)
                            argumentData.values(dim, 0, point) +=
                                    basisData.values(dim, fun, point) *
                                    localCoefficients[fun];
            }
            if (basisDeps & DERIVATIVES)
            {
                std::fill(argumentData.derivatives.begin(),
                          argumentData.derivatives.end(), 0.);
                for (size_t point = 0; point < basisData.derivatives.extent(3); ++point)
                    for (size_t dim = 0; dim < basisData.derivatives.extent(1); ++dim)
                        for (size_t comp = 0; comp < basisData.derivatives.extent(0); ++comp)
                            for (size_t fun = 0; fun < basisData.derivatives.extent(2); ++fun)
                                argumentData.derivatives(comp, dim, 0, point) +=
                                    basisData.derivatives(comp, dim, fun, point) *
                                    localCoefficients[fun];
            }

            // Get geometrical data
            m_rawGeometry->setupGeometry(e, *geometry);
            geometry->getData(trialGeomDeps, localQuadPoints,
                              geomDataPerElement[e]);

            m_trialTransformations->evaluate(argumentData, geomDataPerElement[e],
                                             trialValues);

//            weightedTrialTransfValuesPerElement[e].set_size(transformationCount);
//            for (int transf = 0; transf < transformationCount; ++transf)
//            {
//                const size_t dimCount = trialValues[transf].extent(0);
//                const size_t quadPointCount = trialValues[transf].extent(2);
//                weightedTrialTransfValuesPerElement[e][transf].set_size(
//                            dimCount, quadPointCount);
//                for (size_t point = 0; point < quadPointCount; ++point)
//                    for (size_t dim = 0; dim < dimCount; ++dim)
//                        weightedTrialTransfValuesPerElement[e][transf](dim, point) =
//                                trialValues[transf](dim, 0, point) *
//                                geomDataPerElement[e].integrationElements(point) *
//                                quadWeights[point];
//            } // end of loop over transformations

            const size_t localQuadPointCount = quadWeights.size();

            trialTransfValuesPerElement[e].set_size(transformationCount);
            for (int transf = 0; transf < transformationCount; ++transf)
            {
                const size_t dimCount = trialValues[transf].extent(0);
                assert(trialValues[transf].extent(2) == localQuadPointCount);
                trialTransfValuesPerElement[e][transf].set_size(
                            dimCount, localQuadPointCount);
                for (size_t point = 0; point < localQuadPointCount; ++point)
                    for (size_t dim = 0; dim < dimCount; ++dim)
                        trialTransfValuesPerElement[e][transf](dim, point) =
                                trialValues[transf](dim, 0, point);
            } // end of loop over transformations

            weightsPerElement[e].resize(localQuadPointCount);
            for (size_t point = 0; point < localQuadPointCount; ++point)
                weightsPerElement[e][point] = quadWeights[point] *
                        geomDataPerElement[e].integrationElements(point);

            quadPointCount += quadWeights.size();
        } // end of loop over elements
    } // end of loop over unique bases

    // In the following, weightedTrialExprValuesPerElement[e][transf].extent(1) is used
    // repeatedly as the number of quadrature points in e'th element

    // Now convert std::vectors of arrays into unique big arrays
    // and store them in trialGeomData and weightedTrialTransfValues

    // Fill member matrices of trialGeomData
    if (kernelTrialGeomDeps & GLOBALS)
        trialGeomData.globals.set_size(worldDim, quadPointCount);
    if (kernelTrialGeomDeps & INTEGRATION_ELEMENTS)
        trialGeomData.integrationElements.set_size(quadPointCount);
    if (kernelTrialGeomDeps & NORMALS)
        trialGeomData.normals.set_size(worldDim, quadPointCount);
    if (kernelTrialGeomDeps & JACOBIANS_TRANSPOSED)
        trialGeomData.jacobiansTransposed.set_size(gridDim, worldDim, quadPointCount);
    if (kernelTrialGeomDeps & JACOBIAN_INVERSES_TRANSPOSED)
        trialGeomData.jacobianInversesTransposed.set_size(worldDim, gridDim, quadPointCount);
//    weightedTrialTransfValues.set_size(transformationCount);
//    for (int transf = 0; transf < transformationCount; ++transf)
//        weightedTrialTransfValues[transf].set_size(
//                    m_trialTransformations->resultDimension(transf), quadPointCount);
    trialTransfValues.set_size(transformationCount);
    for (int transf = 0; transf < transformationCount; ++transf)
        trialTransfValues[transf].set_size(
                    m_trialTransformations->resultDimension(transf), quadPointCount);
    weights.resize(quadPointCount);

    for (int e = 0, startCol = 0;
         e < elementCount;
         startCol += trialTransfValuesPerElement[e][0].extent(1), ++e)
    {
        int endCol = startCol + trialTransfValuesPerElement[e][0].extent(1) - 1;
        if (kernelTrialGeomDeps & GLOBALS)
            trialGeomData.globals.cols(startCol, endCol) =
                    geomDataPerElement[e].globals;
        if (kernelTrialGeomDeps & INTEGRATION_ELEMENTS)
            trialGeomData.integrationElements.cols(startCol, endCol) =
                    geomDataPerElement[e].integrationElements;
        if (kernelTrialGeomDeps & NORMALS)
            trialGeomData.normals.cols(startCol, endCol) =
                    geomDataPerElement[e].normals;
        if (kernelTrialGeomDeps & JACOBIANS_TRANSPOSED)
            trialGeomData.jacobiansTransposed.slices(startCol, endCol) =
                    geomDataPerElement[e].jacobiansTransposed;
        if (kernelTrialGeomDeps & JACOBIAN_INVERSES_TRANSPOSED)
            trialGeomData.jacobianInversesTransposed.slices(startCol, endCol) =
                    geomDataPerElement[e].jacobianInversesTransposed;
        for (int transf = 0; transf < transformationCount; ++transf)
            for (size_t point = 0; point < trialTransfValuesPerElement[e][transf].extent(1); ++point)
                for (size_t dim = 0; dim < trialTransfValuesPerElement[e][transf].extent(0); ++dim)
                    trialTransfValues[transf](dim, startCol + point) =
                            trialTransfValuesPerElement[e][transf](dim, point);
        for (size_t point = 0; point < trialTransfValuesPerElement[e][0].extent(1); ++point)
                weights[startCol + point] = weightsPerElement[e][point];
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::quadOrder(
        const Basis<BasisFunctionType>& basis, Region region) const
{
    if (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD)
        return nearFieldQuadOrder(basis);
    else
        return farFieldQuadOrder(basis);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::farFieldQuadOrder(
        const Basis<BasisFunctionType>& basis) const
{
    int elementOrder = (basis.order());
    // Order required for exact quadrature on affine elements with kernel
    // approximated by a polynomial of order identical with that of the basis
    int defaultQuadratureOrder = 2 * elementOrder;
    return m_quadratureOptions.quadratureOrder(defaultQuadratureOrder);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int DefaultEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::nearFieldQuadOrder(
        const Basis<BasisFunctionType>& basis) const
{
    // quick and dirty
    return 2 * farFieldQuadOrder(basis);
}

} // namespace Fiber
