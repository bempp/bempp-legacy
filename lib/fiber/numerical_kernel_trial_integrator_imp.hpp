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

#include "numerical_kernel_trial_integrator.hpp" // To keep IDEs happy

#include "_2d_array.hpp"
#include "_3d_array.hpp"
#include "_4d_array.hpp"

#include "shapeset.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "collection_of_shapeset_transformations.hpp"
#include "geometrical_data.hpp"
#include "collection_of_kernels.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "kernel_trial_integral.hpp"
#include "types.hpp"

#include <cassert>
#include <memory>

namespace Fiber
{

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
NumericalKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
NumericalKernelTrialIntegrator(
        const arma::Mat<CoordinateType>& localQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const arma::Mat<CoordinateType>& points,
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const CollectionOfKernels<KernelType>& kernels,
        const CollectionOfShapesetTransformations<CoordinateType>& trialTransformations,
        const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>& integral) :
    m_localQuadPoints(localQuadPoints),
    m_quadWeights(quadWeights),
    m_points(points),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_kernels(kernels),
    m_trialTransformations(trialTransformations),
    m_integral(integral)
{
    if (localQuadPoints.n_cols != quadWeights.size())
        throw std::invalid_argument("NumericalKernelTrialIntegrator::"
                                    "NumericalKernelTrialIntegrator(): "
                                    "numbers of test points and weights "
                                    "do not match");
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void NumericalKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
integrate(const std::vector<int>& pointIndices,
          int trialElementIndex,
          const Shapeset<BasisFunctionType>& trialShapeset,
          LocalDofIndex localTrialDofIndex,
          const std::vector<arma::Mat<ResultType>*>& result) const
{
    const int quadPointCount = m_localQuadPoints.n_cols;
    const int pointCount = pointIndices.size();
    const int componentCount = m_integral.resultDimension();
    const int trialDofCount = localTrialDofIndex == ALL_DOFS ? trialShapeset.size() : 1;

    if (result.size() != pointCount)
        throw std::invalid_argument(
        "NumericalKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'pointCount' must have the same number "
        "of elements");
    if (pointCount == 0 || quadPointCount == 0)
        return;
    // TODO: in the (pathological) case that quadPointCount == 0 but
    // geometryCount != 0, set elements of result to 0.

    BasisData<BasisFunctionType> trialBasisData;
    GeometricalData<CoordinateType> pointGeomData, trialGeomData;

    size_t trialBasisDeps = 0;
    size_t pointGeomDeps = 0, trialGeomDeps = 0;

    m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernels.addGeometricalDependencies(pointGeomDeps, trialGeomDeps);
    m_integral.addGeometricalDependencies(trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::unique_ptr<Geometry> trialGeometry = m_geometryFactory.make();

    CollectionOf3dArrays<BasisFunctionType> trialValues;
    CollectionOf4dArrays<KernelType> kernelValues;

    for (size_t i = 0; i < result.size(); ++i) {
        assert(result[i]);
        result[i]->set_size(componentCount, trialDofCount);
    }

    m_rawGeometry.setupGeometry(trialElementIndex, *trialGeometry);
    trialShapeset.evaluate(trialBasisDeps, m_localQuadPoints,
                        localTrialDofIndex, trialBasisData);
    trialGeometry->getData(trialGeomDeps, m_localQuadPoints, trialGeomData);
    if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData.domainIndex = m_rawGeometry.domainIndex(trialElementIndex);
    m_trialTransformations.evaluate(trialBasisData, trialGeomData, trialValues);

    // Iterate over the points
    for (int i = 0; i < pointCount; ++i) {
        pointGeomData.globals = m_points.col(pointIndices[i]);
        m_kernels.evaluateOnGrid(pointGeomData, trialGeomData, kernelValues);
        _3dArray<ResultType> result3dView(componentCount, trialDofCount, 1,
                                          result[i]->memptr(), true /* strict */);
        m_integral.evaluateWithPureWeights(
                    trialGeomData, kernelValues, trialValues,
                    m_quadWeights,
                    result3dView);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void NumericalKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
integrate(int pointIndex,
          int componentIndex,
          const std::vector<int>& trialElementIndices,
          const Shapeset<BasisFunctionType>& trialShapeset,
          const std::vector<arma::Mat<ResultType>*>& result) const
{
    const int quadPointCount = m_localQuadPoints.n_cols;
    const int trialElementCount = trialElementIndices.size();
    const int componentCount =
            componentIndex == ALL_COMPONENTS ? m_integral.resultDimension() : 1;
    const int trialDofCount = trialShapeset.size();

    if (result.size() != trialElementCount)
        throw std::invalid_argument(
        "NumericalKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'pointCount' must have the same number "
        "of elements");
    if (trialElementCount == 0 || quadPointCount == 0)
        return;
    // TODO: in the (pathological) case that quadPointCount == 0 but
    // geometryCount != 0, set elements of result to 0.

    BasisData<BasisFunctionType> trialBasisData;
    GeometricalData<CoordinateType> pointGeomData, trialGeomData;

    size_t trialBasisDeps = 0;
    size_t pointGeomDeps = 0, trialGeomDeps = 0;

    m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernels.addGeometricalDependencies(pointGeomDeps, trialGeomDeps);
    m_integral.addGeometricalDependencies(trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::unique_ptr<Geometry> trialGeometry = m_geometryFactory.make();

    CollectionOf3dArrays<BasisFunctionType> trialValues;
    CollectionOf4dArrays<KernelType> kernelValues;

    for (size_t i = 0; i < result.size(); ++i) {
        assert(result[i]);
        result[i]->set_size(componentCount, trialDofCount);
    }
    _3dArray<ResultType> result3d(componentCount, trialDofCount, 1);

    pointGeomData.globals = m_points.col(pointIndex);

    // Iterate over the trial elements
    for (int i = 0; i < trialElementCount; ++i) {
        const int trialElementIndex = trialElementIndices[i];
        m_rawGeometry.setupGeometry(trialElementIndex, *trialGeometry);
        trialShapeset.evaluate(trialBasisDeps, m_localQuadPoints,
                            ALL_DOFS, trialBasisData);
        trialGeometry->getData(trialGeomDeps, m_localQuadPoints, trialGeomData);
        if (trialGeomDeps & DOMAIN_INDEX)
            trialGeomData.domainIndex = m_rawGeometry.domainIndex(trialElementIndex);
        m_trialTransformations.evaluate(trialBasisData, trialGeomData,
                                        trialValues);
        m_kernels.evaluateOnGrid(pointGeomData, trialGeomData, kernelValues);
        // Currently we evaluate all the components of the integral.
        // Later we may consider optimizing and evaluating only one
        // if componentIndex != ALL_COMPONENTS.
        m_integral.evaluateWithPureWeights(
                    trialGeomData, kernelValues, trialValues,
                    m_quadWeights,
                    result3d);
        if (componentIndex == ALL_COMPONENTS)
            for (size_t dof = 0; dof < trialDofCount; ++dof)
                for (size_t component = 0; component < componentCount; ++component)
                    (*result[i])(component, dof) = result3d(component, dof, 0);
        else
            for (size_t dof = 0; dof < trialDofCount; ++dof)
                (*result[i])(0, dof) = result3d(componentIndex, dof, 0);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void NumericalKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
integrate(const std::vector<PointElementIndexPair>& pointElementIndexPairs,
          const Shapeset<BasisFunctionType>& trialShapeset,
          const std::vector<arma::Mat<ResultType>*>& result) const
{
    const int quadPointCount = m_localQuadPoints.n_cols;
    const int pairCount = pointElementIndexPairs.size();
    const int componentCount = m_integral.resultDimension();
    const int trialDofCount = trialShapeset.size();

    if (result.size() != pairCount)
        throw std::invalid_argument(
        "NumericalKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'pointCount' must have the same number "
        "of elements");
    if (pairCount == 0 || quadPointCount == 0)
        return;
    // TODO: in the (pathological) case that quadPointCount == 0 but
    // geometryCount != 0, set elements of result to 0.

    BasisData<BasisFunctionType> trialBasisData;
    GeometricalData<CoordinateType> pointGeomData, trialGeomData;

    size_t trialBasisDeps = 0;
    size_t pointGeomDeps = 0, trialGeomDeps = 0;

    m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernels.addGeometricalDependencies(pointGeomDeps, trialGeomDeps);
    m_integral.addGeometricalDependencies(trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::unique_ptr<Geometry> trialGeometry = m_geometryFactory.make();

    CollectionOf3dArrays<BasisFunctionType> trialValues;
    CollectionOf4dArrays<KernelType> kernelValues;

    for (size_t i = 0; i < result.size(); ++i) {
        assert(result[i]);
        result[i]->set_size(componentCount, trialDofCount);
    }

    // Iterate over the (point, trial element) pairs
    for (int i = 0; i < pairCount; ++i) {
        const int activePointIndex = pointElementIndexPairs[i].first;
        const int activeTrialElementIndex = pointElementIndexPairs[i].second;

        pointGeomData.globals = m_points.col(activePointIndex);
        m_rawGeometry.setupGeometry(activeTrialElementIndex, *trialGeometry);
        trialShapeset.evaluate(trialBasisDeps, m_localQuadPoints,
                            ALL_DOFS, trialBasisData);
        trialGeometry->getData(trialGeomDeps, m_localQuadPoints, trialGeomData);
        if (trialGeomDeps & DOMAIN_INDEX)
            trialGeomData.domainIndex =
                    m_rawGeometry.domainIndex(activeTrialElementIndex);
        m_trialTransformations.evaluate(trialBasisData, trialGeomData, trialValues);
        m_kernels.evaluateOnGrid(pointGeomData, trialGeomData, kernelValues);
        _3dArray<ResultType> result3dView(componentCount, trialDofCount, 1,
                                          result[i]->memptr(), true /* strict */);
        m_integral.evaluateWithPureWeights(
                    trialGeomData, kernelValues, trialValues,
                    m_quadWeights,
                    result3dView);
    }
}

} // namespace Fiber
