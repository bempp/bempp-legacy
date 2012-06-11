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

#include "grid_function.hpp"

#include "assembly_options.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/stl_io.hpp"
#include "../fiber/collection_of_3d_arrays.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/collection_of_basis_transformations.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/function.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/local_assembler_for_grid_functions.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/entity.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer_helper.hpp"
#include "../space/space.hpp"

#include <set>

// TODO: rewrite the constructor of OpenClHandler.
// It should take a bool useOpenCl and *in addition to that* openClOptions.
// The role of the latter should be to e.g. select the device to use
// and other configurable execution parameters.
// If there are no such parameters, OpenClOptions should just be removed.

// Justification: right now there can be a conflict: the user can invoke
// AssemblyOptions::switchToOpenCl() and pass to it an instance of OpenClOptions
// with useOpenCl set to false. This makes no sense.

namespace Bempp
{

namespace
{

template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> reallyCalculateProjections(
        const Space<BasisFunctionType>& space,
        Fiber::LocalAssemblerForGridFunctions<ResultType>& assembler,
        const AssemblyOptions& options)
{
    // TODO: parallelise using TBB (the parameter options will then start be used)

    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = space.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        space.globalDofs(element, testGlobalDofs[elementIndex]);
        it->next();
    }

    // Make a vector of all element indices
    std::vector<int> testIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        testIndices[i] = i;

    // Create the weak form's column vector
    arma::Col<ResultType> result(space.globalDofCount());
    result.fill(0.);

    std::vector<arma::Col<ResultType> > localResult;
    // Evaluate local weak forms
    assembler.evaluateLocalWeakForms(testIndices, localResult);

    // Loop over test indices
    for (int testIndex = 0; testIndex < elementCount; ++testIndex)
        // Add the integrals to appropriate entries in the global weak form
        for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
            result(testGlobalDofs[testIndex][testDof]) +=
                    localResult[testIndex](testDof);

    // Return the vector of projections <phi_i, f>
    return result;
}

/** \brief Calculate projections of the function on test functions from
  the given space. */
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> calculateProjections(
        const Fiber::Function<ResultType>& globalFunction,
        const Space<BasisFunctionType>& space,
        const typename GridFunction<BasisFunctionType, ResultType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& options)
{
    if (!space.dofsAssigned())
        throw std::runtime_error(
                "GridFunction::calculateProjections(): "
                "degrees of freedom of the provided space must be assigned "
                "before calling calculateProjections()");

    // Prepare local assembler
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
    typedef LocalAssemblerConstructionHelper Helper;

    shared_ptr<RawGridGeometry> rawGeometry;
    shared_ptr<GeometryFactory> geometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> testBases;

    Helper::collectGridData(space.grid(),
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelisationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectBases(space, testBases);

    // Get reference to the test basis transformation
    const Fiber::CollectionOfBasisTransformations<CoordinateType>&
            testTransformations = space.shapeFunctionValue();

    typedef Fiber::LocalAssemblerForGridFunctions<ResultType> LocalAssembler;
    std::auto_ptr<LocalAssembler> assembler =
            factory.makeAssemblerForGridFunctions(
                geometryFactory, rawGeometry,
                testBases,
                make_shared_from_ref(testTransformations),
                make_shared_from_ref(globalFunction),
                openClHandler);

    return reallyCalculateProjections(space, *assembler, options);
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromFiberFunction(
        const Space<BasisFunctionType>& space,
        const Fiber::Function<ResultType>& function,
        const typename GridFunction<BasisFunctionType, ResultType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& assemblyOptions)
{
    arma::Col<ResultType> projections =
            calculateProjections(function, space, factory, assemblyOptions);
    return gridFunctionFromProjections(space, projections);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromCoefficients(const Space<BasisFunctionType>& space,
                             const arma::Col<ResultType>& coefficients)
{
    typedef GridFunction<BasisFunctionType, ResultType> GF;

    if (!space.dofsAssigned())
        throw std::runtime_error(
                "GridFunction::GridFunction(): "
                "degrees of freedom of the provided space must be assigned "
                "beforehand");
    if (space.globalDofCount() != coefficients.n_rows)
        throw std::runtime_error(
                "GridFunction::GridFunction(): "
                "dimension of coefficients does not match the number of global "
                "DOFs in the provided function space");

    arma::Col<ResultType> projections;
    space.applyMassMatrix(coefficients, projections);
    return GF(space, coefficients, projections);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
gridFunctionFromProjections(const Space<BasisFunctionType>& space,
                            const arma::Col<ResultType>& projections)
{
    typedef GridFunction<BasisFunctionType, ResultType> GF;

    if (!space.dofsAssigned())
        throw std::runtime_error(
                "GridFunction::GridFunction(): "
                "degrees of freedom of the provided space must be assigned "
                "beforehand");
    if (space.globalDofCount() != projections.n_rows)
        throw std::runtime_error(
                "GridFunction::GridFunction(): "
                "dimension of projections does not match the number of global "
                "DOFs in the provided function space");

    arma::Col<ResultType> coefficients;
    space.applyInverseMassMatrix(projections, coefficients);
    return GF(space, coefficients, projections);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const Space<BasisFunctionType>& space,
        const arma::Col<ResultType>& coefficients,
        const arma::Col<ResultType>& projections) :
    m_space(space), m_coefficients(coefficients), m_projections(projections)
{
    // TODO: assert that dimensions are correct
}

template <typename BasisFunctionType, typename ResultType>
const Grid& GridFunction<BasisFunctionType, ResultType>::grid() const
{
    return m_space.grid();
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>& GridFunction<BasisFunctionType, ResultType>::space() const
{
    return m_space;
}

template <typename BasisFunctionType, typename ResultType>
int GridFunction<BasisFunctionType, ResultType>::codomainDimension() const
{
    return m_space.codomainDimension();
}

template <typename BasisFunctionType, typename ResultType>
const arma::Col<ResultType>& GridFunction<BasisFunctionType, ResultType>::coefficients() const
{
    return m_coefficients;
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setCoefficients(
        const arma::Col<ResultType>& coeffs)
{
    if (coeffs.n_rows != m_space.globalDofCount())
        throw std::invalid_argument(
                "GridFunction::setCoefficients(): dimension of the provided "
                "vector does not match the number of global DOFs");
    m_coefficients = coeffs;
    m_space.applyMassMatrix(coeffs, m_projections);
}

template <typename BasisFunctionType, typename ResultType>
const arma::Col<ResultType>& GridFunction<BasisFunctionType, ResultType>::projections() const
{
    return m_projections;
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setProjections(
        const arma::Col<ResultType>& projects)
{
    if (projects.n_rows != m_space.globalDofCount())
        throw std::invalid_argument(
                "GridFunction::setCoefficients(): dimension of the provided "
                "vector does not match the number of global DOFs");
    m_projections = projects;
    m_space.applyInverseMassMatrix(projects, m_coefficients);
}

// Redundant, in fact -- can be obtained directly from Space
template <typename BasisFunctionType, typename ResultType>
const Fiber::Basis<BasisFunctionType>& GridFunction<BasisFunctionType, ResultType>::basis(
        const Entity<0>& element) const
{
    return m_space.basis(element);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::getLocalCoefficients(
        const Entity<0>& element, std::vector<ResultType>& coeffs) const
{
    std::vector<GlobalDofIndex> gdofIndices;
    m_space.globalDofs(element, gdofIndices);
    const int gdofCount = gdofIndices.size();
    coeffs.resize(gdofCount);
    for (int i = 0; i < gdofCount; ++i)
        coeffs[i] = m_coefficients(gdofIndices[i]);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::exportToVtk(
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType) const
{
    arma::Mat<ResultType> data;
    evaluateAtSpecialPoints(dataType, data);

    std::auto_ptr<GridView> view = m_space.grid().leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();

    exportSingleDataSetToVtk(*vtkWriter, data, dataType, dataLabel,
                             fileNamesBase, filesPath, outputType);
}


template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::evaluateAtSpecialPoints(
        VtkWriter::DataType dataType, arma::Mat<ResultType>& result) const
{
    if (dataType != VtkWriter::CELL_DATA && dataType != VtkWriter::VERTEX_DATA)
        throw std::invalid_argument("GridFunction::evaluateAtSpecialPoints(): "
                                    "invalid data type");

    const Grid& grid = m_space.grid();
    const int gridDim = grid.dim();
    const int elementCodim = 0;
    const int vertexCodim = grid.dim();
    const int codomainDim = codomainDimension();

    std::auto_ptr<GridView> view = grid.leafView();
    const int elementCount = view->entityCount(elementCodim);
    const int vertexCount = view->entityCount(vertexCodim);

    result.set_size(codomainDimension(),
                    dataType == VtkWriter::CELL_DATA ? elementCount : vertexCount);
    result.fill(0.);

    // Number of elements contributing to each column in result
    // (this will be greater than 1 for VERTEX_DATA)
    std::vector<int> multiplicities(vertexCount);
    std::fill(multiplicities.begin(), multiplicities.end(), 0);

    // Gather geometric data
    Fiber::RawGridGeometry<CoordinateType> rawGeometry(grid.dim(), grid.dimWorld());
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            grid.elementGeometryFactory();
    std::auto_ptr<typename GeometryFactory::Geometry> geometry(
                geometryFactory->make());
    Fiber::GeometricalData<CoordinateType> geomData;

    // For each element, get its basis and corner count (this is sufficient
    // to identify its geometry) as well as its local coefficients
    typedef std::pair<const Fiber::Basis<BasisFunctionType>*, int> BasisAndCornerCount;
    typedef std::vector<BasisAndCornerCount> BasisAndCornerCountVector;
    BasisAndCornerCountVector basesAndCornerCounts(elementCount);
    std::vector<std::vector<ResultType> > localCoefficients(elementCount);
    {
        std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
        for (int e = 0; e < elementCount; ++e) {
            const Entity<0>& element = it->entity();
            basesAndCornerCounts[e] = BasisAndCornerCount(
                        &m_space.basis(element), rawGeometry.elementCornerCount(e));
            getLocalCoefficients(element, localCoefficients[e]);
            it->next();
        }
    }

    typedef std::set<BasisAndCornerCount> BasisAndCornerCountSet;
    BasisAndCornerCountSet uniqueBasesAndCornerCounts(
                basesAndCornerCounts.begin(), basesAndCornerCounts.end());

    // Find out which basis data need to be calculated
    int basisDeps = 0, geomDeps = 0;
    // Find out which geometrical data need to be calculated, in addition
    // to those needed by the kernel
    const Fiber::CollectionOfBasisTransformations<CoordinateType>& transformations =
            m_space.shapeFunctionValue();
    assert(codomainDim == transformations.resultDimension(0));
    transformations.addDependencies(basisDeps, geomDeps);

    // Loop over unique combinations of basis and element corner count
    typedef typename BasisAndCornerCountSet::const_iterator
            BasisAndCornerCountSetConstIt;
    for (BasisAndCornerCountSetConstIt it = uniqueBasesAndCornerCounts.begin();
         it != uniqueBasesAndCornerCounts.end(); ++it) {
        const BasisAndCornerCount& activeBasisAndCornerCount = *it;
        const Fiber::Basis<BasisFunctionType>& activeBasis =
                *activeBasisAndCornerCount.first;
        int activeCornerCount = activeBasisAndCornerCount.second;

        // Set the local coordinates of either all vertices or the barycentre
        // of the active element type
        arma::Mat<CoordinateType> local;
        if (dataType == VtkWriter::CELL_DATA) {
            local.set_size(gridDim, 1);

            // We could actually use Dune for these assignements
            if (gridDim == 1 && activeCornerCount == 2) {
                // linear segment
                local(0, 0) = 0.5;
            } else if (gridDim == 2 && activeCornerCount == 3) {
                // triangle
                local(0, 0) = 1./3.;
                local(1, 0) = 1./3.;
            } else if (gridDim == 2 && activeCornerCount == 4) {
                // quadrilateral
                local(0, 0) = 0.5;
                local(1, 0) = 0.5;
            } else
                throw std::runtime_error("GridFunction::evaluateAtVertices(): "
                                         "unsupported element type");
        } else { // VERTEX_DATA
            local.set_size(gridDim, activeCornerCount);

            // We could actually use Dune for these assignements
            if (gridDim == 1 && activeCornerCount == 2) {
                // linear segment
                local(0, 0) = 0.;
                local(0, 1) = 1.;
            } else if (gridDim == 2 && activeCornerCount == 3) {
                // triangle
                local.fill(0.);
                local(0, 1) = 1.;
                local(1, 2) = 1.;
            } else if (gridDim == 2 && activeCornerCount == 4) {
                // quadrilateral
                local.fill(0.);
                local(0, 1) = 1.;
                local(1, 2) = 1.;
                local(0, 3) = 1.;
                local(1, 3) = 1.;
            } else
                throw std::runtime_error("GridFunction::evaluateAtVertices(): "
                                         "unsupported element type");
        }

        // Get basis data
        Fiber::BasisData<BasisFunctionType> basisData;
        activeBasis.evaluate(basisDeps, local, ALL_DOFS, basisData);

        Fiber::BasisData<ResultType> functionData;
        if (basisDeps & Fiber::VALUES)
            functionData.values.set_size(basisData.values.n_rows,
                                         1, // just one function
                                         basisData.values.n_slices);
        if (basisDeps & Fiber::DERIVATIVES)
            functionData.derivatives.set_size(basisData.derivatives.extent(0),
                                              basisData.derivatives.extent(1),
                                              1, // just one function
                                              basisData.derivatives.extent(3));
        Fiber::CollectionOf3dArrays<ResultType> functionValues;

        // Loop over elements and process those that use the active basis
        for (int e = 0; e < elementCount; ++e) {
            if (basesAndCornerCounts[e].first != &activeBasis)
                continue;

            // Local coefficients of the argument in the current element
            const std::vector<ResultType>& activeLocalCoefficients =
                    localCoefficients[e];

            // Calculate the function's values and/or derivatives
            // at the requested points in the current element
            if (basisDeps & Fiber::VALUES) {
                functionData.values.fill(0.);
                for (int point = 0; point < basisData.values.n_slices; ++point)
                    for (int dim = 0; dim < basisData.values.n_rows; ++dim)
                        for (int fun = 0; fun < basisData.values.n_cols; ++fun)
                            functionData.values(dim, 0, point) +=
                                    basisData.values(dim, fun, point) *
                                    activeLocalCoefficients[fun];
            }
            if (basisDeps & Fiber::DERIVATIVES) {
                std::fill(functionData.derivatives.begin(),
                          functionData.derivatives.end(), 0.);
                for (int point = 0; point < basisData.derivatives.extent(3); ++point)
                    for (int dim = 0; dim < basisData.derivatives.extent(1); ++dim)
                        for (int comp = 0; comp < basisData.derivatives.extent(0); ++comp)
                            for (int fun = 0; fun < basisData.derivatives.extent(2); ++fun)
                                functionData.derivatives(comp, dim, 0, point) +=
                                        basisData.derivatives(comp, dim, fun, point) *
                                        activeLocalCoefficients[fun];
            }

            // Get geometrical data
            rawGeometry.setupGeometry(e, *geometry);
            geometry->getData(geomDeps, local, geomData);

            transformations.evaluate(functionData, geomData, functionValues);
            assert(functionValues[0].extent(1) == 1); // one function

            if (dataType == VtkWriter::CELL_DATA)
                for (int dim = 0; dim < codomainDim; ++dim)
                    result(dim, e) = functionValues[0](     // array index
                                                       dim, // component
                                                       0,   // function index
                                                       0);  // point index
            else { // VERTEX_DATA
                // Add the calculated values to the columns of the result array
                // corresponding to the active element's vertices
                for (int c = 0; c < activeCornerCount; ++c) {
                    int vertexIndex = rawGeometry.elementCornerIndices()(c, e);
                    for (int dim = 0; dim < codomainDim; ++dim)
                        result(dim, vertexIndex) += functionValues[0](dim, 0, c);
                    ++multiplicities[vertexIndex];
                }
            }
        } // end of loop over elements
    } // end of loop over unique combinations of basis and corner count

    // Take average of the vertex values obtained in each of the adjacent elements
    if (dataType == VtkWriter::VERTEX_DATA)
        for (int v = 0; v < vertexCount; ++v)
            result.col(v) /= multiplicities[v];
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2)
{
    if (&g1.space() != &g2.space())
        throw std::runtime_error("GridFunction::operator+(): spaces don't match");
    return GridFunction<BasisFunctionType, ResultType>(
                g1.space(),
                g1.coefficients() + g2.coefficients(),
                g1.projections() + g2.projections());
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2)
{
    if (&g1.space() != &g2.space())
        throw std::runtime_error("GridFunction::operator-(): spaces don't match");
    return GridFunction<BasisFunctionType, ResultType>(
                g1.space(),
                g1.coefficients() - g2.coefficients(),
                g1.projections() - g2.projections());
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar)
{
    return GridFunction<BasisFunctionType, ResultType>(
                g1.space(),
                static_cast<ResultType>(scalar) * g1.coefficients(),
                static_cast<ResultType>(scalar) * g1.projections());
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    GridFunction<BasisFunctionType, ResultType> >::type
operator*(
        const ScalarType& scalar, const GridFunction<BasisFunctionType, ResultType>& g2)
{
    return g2 * scalar;
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator/(
        const GridFunction<BasisFunctionType, ResultType>& g1, const ScalarType& scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("GridFunction::operator/(): Divide by zero");
    return (static_cast<ScalarType>(1.) / scalar) * g1;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template GridFunction<BASIS, RESULT> \
    gridFunctionFromFiberFunction( \
    const Space<BASIS>& space, \
    const Fiber::Function<RESULT>& function, \
    const GridFunction<BASIS, RESULT>::LocalAssemblerFactory& factory, \
    const AssemblyOptions& assemblyOptions); \
    \
    template GridFunction<BASIS, RESULT> \
    gridFunctionFromCoefficients(const Space<BASIS>& space, \
    const arma::Col<RESULT>& coefficients); \
    \
    template GridFunction<BASIS, RESULT> \
    gridFunctionFromProjections(const Space<BASIS>& space, \
    const arma::Col<RESULT>& coefficients); \
    \
    template GridFunction<BASIS, RESULT> operator+( \
    const GridFunction<BASIS, RESULT>& op1, \
    const GridFunction<BASIS, RESULT>& op2); \
    \
    template GridFunction<BASIS, RESULT> operator-( \
    const GridFunction<BASIS, RESULT>& op1, \
    const GridFunction<BASIS, RESULT>& op2)
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(BASIS, RESULT, SCALAR) \
    template GridFunction<BASIS, RESULT> operator*( \
    const GridFunction<BASIS, RESULT>& op, const SCALAR& scalar); \
    template GridFunction<BASIS, RESULT> operator*( \
    const SCALAR& scalar, const GridFunction<BASIS, RESULT>& op); \
    template GridFunction<BASIS, RESULT> operator/( \
    const GridFunction<BASIS, RESULT>& op, const SCALAR& scalar)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(
        float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, float, double);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && (defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(
        float, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, std::complex<double>);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(
        std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(
        double, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && (defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(
        double, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(
        std::complex<double>, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, std::complex<double>);
#endif

} // namespace Bempp
