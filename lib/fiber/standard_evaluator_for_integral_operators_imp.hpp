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

#include "standard_evaluator_for_integral_operators.hpp" // keep IDEs happy

#include "basis.hpp"
#include "basis_data.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "kernel.hpp"
#include "numerical_quadrature.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"

namespace Fiber
{

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::StandardEvaluatorForIntegralOperators(
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const std::vector<const Basis<BasisFunctionType>*>& trialBases,
        const Kernel<KernelType>& kernel,
        const Expression<CoordinateType>& trialExpression,
        const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
        const OpenClHandler<CoordinateType, int>& openClHandler,
        const QuadratureOptions& quadratureOptions) :
    m_geometryFactory(geometryFactory), m_rawGeometry(rawGeometry),
    m_trialBases(trialBases), m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_argumentLocalCoefficients(argumentLocalCoefficients),
    m_openClHandler(openClHandler), m_quadratureOptions(quadratureOptions)
{
    const int elementCount = rawGeometry.elementCount();
    if (!rawGeometry.auxData().is_empty() &&
            rawGeometry.auxData().n_cols != elementCount)
        throw std::invalid_argument(
                "StandardEvaluatorForIntegralOperators::"
                "StandardEvaluatorForIntegralOperators(): "
                "number of columns of auxData must match that of "
                "elementCornerIndices");
    if (trialBases.size() != elementCount)
        throw std::invalid_argument(
                "StandardEvaluatorForIntegralOperators::"
                "StandardEvaluatorForIntegralOperators(): "
                "size of testBases must match the number of columns of "
                "elementCornerIndices");

    // Assert that the kernel tensor dimensions are compatible
    // with the number of components of the argument

    const int kernelRowCount = m_kernel.codomainDimension();
    const int kernelColCount = m_kernel.domainDimension();
    const int argumentComponentCount = m_trialExpression.codomainDimension();

    const bool scalarKernel = (kernelRowCount == 1 && kernelColCount == 1);
    if (!scalarKernel)
        if (kernelColCount != argumentComponentCount)
            throw std::invalid_argument(
                    "StandardEvaluatorForIntegralOperators::"
                    "StandardEvaluatorForIntegralOperators(): "
                    "incompatible dimensions of kernel and argument");

    // Cache "trial data" such as the values of the argument at quadrature
    // points, vectors normal to surface at quadrature points etc.
    cacheTrialData();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::evaluate(
        Region region,
        const arma::Mat<CoordinateType>& points, arma::Mat<ResultType>& result) const
{
    const int pointCount = points.n_cols;
    const int worldDim = points.n_rows;
    if (worldDim != m_kernel.worldDimension())
        throw std::invalid_argument(
                "StandardEvaluatorForIntegralOperators::evaluateFarField(): "
                "points have incorrect number of components");

    const int kernelRowCount = m_kernel.codomainDimension();
    const int kernelColCount = m_kernel.domainDimension();
    const int argumentComponentCount = m_trialExpression.codomainDimension();

    const bool scalarKernel = (kernelRowCount == 1 && kernelColCount == 1);
    const int outputComponentCount =
            scalarKernel ? argumentComponentCount : kernelRowCount;

    result.set_size(outputComponentCount, pointCount);
    result.fill(0.);

    const Fiber::GeometricalData<CoordinateType>& trialGeomData =
            (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD) ?
                m_nearFieldTrialGeomData :
                m_farFieldTrialGeomData;
    const arma::Mat<ResultType>& weightedTrialExprValues =
            (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD) ?
                m_nearFieldWeightedTrialExprValues :
                m_farFieldWeightedTrialExprValues;

    const int quadPointCount = weightedTrialExprValues.n_cols;

    // Do things in chunks of 96 points -- in order to avoid creating
    // too large arrays of kernel values
    const int chunkSize = 96;
    Fiber::Array4D<KernelType> kernelValues;
    Fiber::GeometricalData<CoordinateType> testGeomData;
    for (int start = 0; start < pointCount; start += chunkSize)
    {
        int end = std::min(start + chunkSize - 1, pointCount - 1);
        testGeomData.globals = points.cols(start, end);
        m_kernel.evaluateOnGrid(testGeomData, trialGeomData, kernelValues);
        if (scalarKernel)
            for (int evalPoint = start; evalPoint <= end; ++evalPoint) // "<=" intended!
                for (int dim = 0; dim < outputComponentCount; ++dim)
                    for (int quadPoint = 0; quadPoint < quadPointCount; ++quadPoint)
                        result(dim, evalPoint) +=
                                kernelValues(0, evalPoint - start, 0, quadPoint) *
                                weightedTrialExprValues(dim, quadPoint);
        else
            for (int evalPoint = start; evalPoint <= end; ++evalPoint) // "<=" intended!
                for (int kernelRow = 0; kernelRow < kernelRowCount; ++kernelRow)
                    for (int kernelCol = 0; kernelCol < kernelColCount; ++kernelCol)
                        for (int quadPoint = 0; quadPoint < quadPointCount; ++quadPoint)
                            result(kernelRow, evalPoint) +=
                                    kernelValues(kernelRow, evalPoint - start, kernelCol, quadPoint) *
                                    weightedTrialExprValues(kernelCol, quadPoint);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::cacheTrialData()
{
    int testGeomDeps = 0, trialGeomDeps = 0;
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
    if (testGeomDeps != 0 && testGeomDeps != Fiber::GLOBALS)
        throw std::runtime_error(
                "StandardEvaluatorForIntegralOperators::cacheTrialData(): "
                "integral operators whose kernels depend on other test data "
                "than global coordinates cannot be evaluated off-surface");

    calcTrialData(EvaluatorForIntegralOperators<ResultType>::FAR_FIELD,
                  trialGeomDeps, m_farFieldTrialGeomData,
                  m_farFieldWeightedTrialExprValues);
    // near field currently not used
    calcTrialData(EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD,
                  trialGeomDeps, m_nearFieldTrialGeomData,
                  m_nearFieldWeightedTrialExprValues);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::calcTrialData(
        Region region,
        int kernelTrialGeomDeps,
        Fiber::GeometricalData<CoordinateType>& trialGeomData,
        arma::Mat<ResultType>& weightedTrialExprValues) const
{
    const int elementCount = m_rawGeometry.elementCount();
    const int worldDim = m_rawGeometry.worldDimension();
    const int gridDim = m_rawGeometry.gridDimension();
    const int trialExprComponentCount = m_trialExpression.codomainDimension();

    // Find out which basis data need to be calculated
    int basisDeps = 0;
    // Find out which geometrical data need to be calculated, in addition
    // to those needed by the kernel
    int trialGeomDeps = kernelTrialGeomDeps;
    m_trialExpression.addDependencies(basisDeps, trialGeomDeps);
    trialGeomDeps |= INTEGRATION_ELEMENTS;

    // Initialise the geometry factory
    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory.make());

    // Find all unique trial bases
    // Set of unique quadrature variants
    typedef std::set<const Basis<BasisFunctionType>*> BasisSet;
    BasisSet uniqueTrialBases(m_trialBases.begin(), m_trialBases.end());

    // Initialise temporary (per-element) data containers
    std::vector<GeometricalData<CoordinateType> > geomDataPerElement(elementCount);
    std::vector<arma::Mat<ResultType> > weightedTrialExprValuesPerElement(elementCount);

    for (typename BasisSet::const_iterator it = uniqueTrialBases.begin();
         it != uniqueTrialBases.end(); ++it)
    {
        const Basis<BasisFunctionType>& activeBasis = *(*it);
        int order = quadOrder(activeBasis, region);

        // Find out the element type
        int elementCornerCount = -1;
        for (int e = 0; e < elementCount; ++e)
            if (m_trialBases[e] == &activeBasis)
            {
                elementCornerCount = m_rawGeometry.elementCornerCount(e);
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
            argumentData.values.set_size(basisData.values.n_rows,
                                         1, // just one function
                                         basisData.values.n_slices);
        if (basisDeps & DERIVATIVES)
            argumentData.derivatives.set_size(basisData.derivatives.extent(0),
                                              basisData.derivatives.extent(1),
                                              1, // just one function
                                              basisData.derivatives.extent(3));

        // Loop over elements and process those that use the active basis
        for (int e = 0; e < elementCount; ++e)
        {
            if (m_trialBases[e] != &activeBasis)
                continue;

            // Local coefficients of the argument in the current element
            const std::vector<ResultType>& localCoefficients =
                    m_argumentLocalCoefficients[e];

            // Calculate the argument function's values and/or derivatives
            // at quadrature points in the current element
            if (basisDeps & VALUES)
            {
                argumentData.values.fill(0.);
                for (int point = 0; point < basisData.values.n_slices; ++point)
                    for (int dim = 0; dim < basisData.values.n_rows; ++dim)
                        for (int fun = 0; fun < basisData.values.n_cols; ++fun)
                            argumentData.values(dim, 0, point) +=
                                    basisData.values(dim, fun, point) *
                                    localCoefficients[fun];
            }
            if (basisDeps & DERIVATIVES)
            {
                std::fill(argumentData.derivatives.begin(),
                          argumentData.derivatives.end(), 0.);
                for (int point = 0; point < basisData.derivatives.extent(3); ++point)
                    for (int dim = 0; dim < basisData.derivatives.extent(1); ++dim)
                        for (int comp = 0; comp < basisData.derivatives.extent(0); ++comp)
                            for (int fun = 0; fun < basisData.derivatives.extent(2); ++fun)
                                argumentData.derivatives(comp, dim, 0, point) +=
                                    basisData.derivatives(comp, dim, fun, point) *
                                    localCoefficients[fun];
            }

            // Get geometrical data
            m_rawGeometry.setupGeometry(e, *geometry);
            geometry->getData(trialGeomDeps, localQuadPoints,
                              geomDataPerElement[e]);

            arma::Cube<ResultType> trialValues;
            m_trialExpression.evaluate(argumentData, geomDataPerElement[e],
                                       trialValues);

            assert(trialValues.n_cols == 1);
            weightedTrialExprValuesPerElement[e].set_size(trialValues.n_rows,
                                                          trialValues.n_slices);
            for (int point = 0; point < trialValues.n_slices; ++point)
                for (int dim = 0; dim < trialValues.n_rows; ++dim)
                    weightedTrialExprValuesPerElement[e](dim, point) =
                            trialValues(dim, 0, point) *
                            geomDataPerElement[e].integrationElements(point) *
                            quadWeights[point];
        } // end of loop over elements
    } // end of loop over unique bases

    // In the following, weightedTrialExprValuesPerElement[e].n_cols is used
    // repeatedly as the number of quadrature points in e'th element

    // Count quadrature points
    int quadPointCount = 0;
    for (int e = 0; e < elementCount; ++e)
        quadPointCount += weightedTrialExprValuesPerElement[e].n_cols;

    // Now convert std::vectors of arrays into unique big arrays
    // and store them in trialGeomData

    // Fill member matrices of trialGeomData
    if (kernelTrialGeomDeps & Fiber::GLOBALS)
        trialGeomData.globals.set_size(worldDim, quadPointCount);
    if (kernelTrialGeomDeps & Fiber::INTEGRATION_ELEMENTS)
        trialGeomData.integrationElements.set_size(quadPointCount);
    if (kernelTrialGeomDeps & Fiber::NORMALS)
        trialGeomData.normals.set_size(worldDim, quadPointCount);
    if (kernelTrialGeomDeps & Fiber::JACOBIANS_TRANSPOSED)
        trialGeomData.jacobiansTransposed.set_size(gridDim, worldDim, quadPointCount);
    if (kernelTrialGeomDeps & Fiber::JACOBIAN_INVERSES_TRANSPOSED)
        trialGeomData.jacobianInversesTransposed.set_size(worldDim, gridDim, quadPointCount);
    weightedTrialExprValues.set_size(trialExprComponentCount, quadPointCount);

    for (int e = 0, startCol = 0;
         e < elementCount;
         ++e, startCol += weightedTrialExprValuesPerElement[e].n_cols)
    {
        int endCol = startCol + weightedTrialExprValuesPerElement[e].n_cols - 1;
        if (kernelTrialGeomDeps & Fiber::GLOBALS)
            trialGeomData.globals.cols(startCol, endCol) =
                    geomDataPerElement[e].globals;
        if (kernelTrialGeomDeps & Fiber::INTEGRATION_ELEMENTS)
            trialGeomData.integrationElements.cols(startCol, endCol) =
                    geomDataPerElement[e].integrationElements;
        if (kernelTrialGeomDeps & Fiber::NORMALS)
            trialGeomData.normals.cols(startCol, endCol) =
                    geomDataPerElement[e].normals;
        if (kernelTrialGeomDeps & Fiber::JACOBIANS_TRANSPOSED)
            trialGeomData.jacobiansTransposed.slices(startCol, endCol) =
                    geomDataPerElement[e].jacobiansTransposed;
        if (kernelTrialGeomDeps & Fiber::JACOBIAN_INVERSES_TRANSPOSED)
            trialGeomData.jacobianInversesTransposed.slices(startCol, endCol) =
                    geomDataPerElement[e].jacobianInversesTransposed;
        weightedTrialExprValues.cols(startCol, endCol) =
                weightedTrialExprValuesPerElement[e];
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::quadOrder(
        const Fiber::Basis<BasisFunctionType>& basis, Region region) const
{
    if (region == EvaluatorForIntegralOperators<ResultType>::NEAR_FIELD)
        return nearFieldQuadOrder(basis);
    else
        return farFieldQuadOrder(basis);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::farFieldQuadOrder(
        const Fiber::Basis<BasisFunctionType>& basis) const
{
    // const QuadratureOptions& options = m_accuracyOptions.farField;

    if (m_quadratureOptions.mode == QuadratureOptions::EXACT_ORDER)
        return m_quadratureOptions.order;
    else
    {
        // Order required for exact quadrature on affine elements with kernel
        // approximated by a polynomial of order identical with that of
        // the basis
        int elementOrder = (basis.order());
        int minimumOrder = ((2 * elementOrder + 1) + 1) / 2;
        return minimumOrder + m_quadratureOptions.orderIncrement;
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int StandardEvaluatorForIntegralOperators<BasisFunctionType, KernelType,
ResultType, GeometryFactory>::nearFieldQuadOrder(
        const Fiber::Basis<BasisFunctionType>& basis) const
{
    // quick and dirty
    return 2 * farFieldQuadOrder(basis);
}

} // namespace Fiber
