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

#include "abstract_boundary_operator_pseudoinverse.hpp"
#include "assembly_options.hpp"
#include "boundary_operator.hpp"
#include "context.hpp"
#include "discrete_boundary_operator.hpp"
#include "identity_operator.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/complex_aux.hpp"
#include "../common/deprecated.hpp"
#include "../fiber/collection_of_3d_arrays.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/collection_of_basis_transformations.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/function.hpp"
#include "../fiber/local_assembler_for_grid_functions.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/entity.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer_helper.hpp"
#include "../space/space.hpp"
#include "identity_operator.hpp"
#include "../io/gmsh.hpp"

#include <boost/array.hpp>
#include <fstream>
#include <set>
#include <sstream>

namespace Bempp
{

// Internal routines

namespace
{

template <typename BasisFunctionType, typename ResultType>
shared_ptr<arma::Col<ResultType> > reallyCalculateProjections(
        const Space<BasisFunctionType>& dualSpace,
        Fiber::LocalAssemblerForGridFunctions<ResultType>& assembler,
        const AssemblyOptions& options)
{
    // TODO: parallelise using TBB (the parameter options will then start be used)

    // Get the grid's leaf view so that we can iterate over elements
    const GridView& view = dualSpace.gridView();
    const size_t elementCount = view.entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs(elementCount);
    std::vector<std::vector<BasisFunctionType> > testLocalDofWeights(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view.elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        dualSpace.getGlobalDofs(element, testGlobalDofs[elementIndex],
            testLocalDofWeights[elementIndex]);
        it->next();
    }

    // Make a vector of all element indices
    std::vector<int> testIndices(elementCount);
    for (size_t i = 0; i < elementCount; ++i)
        testIndices[i] = i;

    // Create the weak form's column vector
    shared_ptr<arma::Col<ResultType> > result(
                new arma::Col<ResultType>(dualSpace.globalDofCount()));
    result->fill(0.);

    std::vector<arma::Col<ResultType> > localResult;
    // Evaluate local weak forms
    assembler.evaluateLocalWeakForms(testIndices, localResult);

    // Loop over test indices
    for (size_t testIndex = 0; testIndex < elementCount; ++testIndex)
        // Add the integrals to appropriate entries in the global weak form
        for (size_t testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof) {
            int testGlobalDof = testGlobalDofs[testIndex][testDof];
            if (testGlobalDof >= 0) // if it's negative, it means that this
                                    // local dof is constrained (not used)
                (*result)(testGlobalDof) +=
                    conj(testLocalDofWeights[testIndex][testDof]) *
                    localResult[testIndex](testDof);
        }

    // Return the vector of projections <phi_i, f>
    return result;
}

/** \brief Calculate projections of the function on the basis functions of
  the given dual space. */
template <typename BasisFunctionType, typename ResultType>
shared_ptr<arma::Col<ResultType> > calculateProjections(
        const Context<BasisFunctionType, ResultType>& context,
        const Function<ResultType>& globalFunction,
        const Space<BasisFunctionType>& dualSpace)
{
    const AssemblyOptions& options = context.assemblyOptions();

    // Prepare local assembler
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Shapeset<BasisFunctionType>*> ShapesetPtrVector;
    typedef LocalAssemblerConstructionHelper Helper;

    shared_ptr<RawGridGeometry> rawGeometry;
    shared_ptr<GeometryFactory> geometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<ShapesetPtrVector> testShapesets;

    Helper::collectGridData(dualSpace,
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectShapesets(dualSpace, testShapesets);

    // Get reference to the test shapeset transformation
    const Fiber::CollectionOfShapesetTransformations<CoordinateType>&
            testTransformations = dualSpace.basisFunctionValue();

    typedef Fiber::LocalAssemblerForGridFunctions<ResultType> LocalAssembler;
    std::auto_ptr<LocalAssembler> assembler =
            context.quadStrategy()->makeAssemblerForGridFunctions(
                geometryFactory, rawGeometry,
                testShapesets,
                make_shared_from_ref(testTransformations),
                make_shared_from_ref(globalFunction),
                openClHandler);

    return reallyCalculateProjections(dualSpace, *assembler, options);
}

/** \brief Evaluate the function at the interpolation points of the chosen space. */
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType> interpolate(
        const Function<ResultType>& globalFunction,
        const Space<BasisFunctionType>& space)
{
    size_t deps = 0;
    globalFunction.addGeometricalDependencies(deps);
    if (deps & ~(Fiber::GLOBALS | Fiber::NORMALS))
        throw std::invalid_argument(
                "interpolate(): functions to be interpolated "
                "must not depend on any geometrical data besides global "
                "coordinates and normal vectors");

    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    Fiber::GeometricalData<CoordinateType> geomData;
    if (deps & Fiber::GLOBALS) {
        space.getGlobalDofInterpolationPoints(geomData.globals);
    }
    if (deps & Fiber::NORMALS) {
        space.getNormalsAtGlobalDofInterpolationPoints(geomData.normals);
    }
    arma::Mat<ResultType> values;
    globalFunction.evaluate(geomData, values);

    const size_t componentCount = values.n_rows;
    const size_t pointCount = values.n_cols;

    arma::Mat<CoordinateType> directions;
    space.getGlobalDofInterpolationDirections(directions);
    assert(directions.n_rows == values.n_rows);
    assert(directions.n_cols == pointCount);

    arma::Col<ResultType> result(pointCount);
    result.fill(0);
    for (size_t p = 0; p < pointCount; ++p)
        for (size_t d = 0; d < componentCount; ++d)
            result(p) += values(d, p) * directions(d, p);
    return result;
}

} // namespace

// Recommended constructors

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction()
{
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const arma::Col<ResultType>& coefficients)
{
    initializeFromCoefficients(context, space, coefficients);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
        const arma::Col<ResultType>& projections)
{

    initializeFromProjections(context, space, dualSpace, projections);
    m_dualSpace = dualSpace;
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
        const Function<ResultType>& function,
        ConstructionMode mode) :
    m_context(context), m_space(space), m_dualSpace(dualSpace)
{
    if (!context)
        throw std::invalid_argument(
                "GridFunction::GridFunction(): context must not be null");
    if (!space)
        throw std::invalid_argument(
                "GridFunction::GridFunction(): space must not be null");
    if (!dualSpace)
        throw std::invalid_argument(
                "GridFunction::GridFunction(): dualSpace must not be null");

    if (space->codomainDimension() != dualSpace->codomainDimension())
        throw std::invalid_argument(
                "GridFunction::GridFunction(): "
                "functions from 'space' and 'dualSpace' have a different "
                "number of components");
    if (function.codomainDimension() != space->codomainDimension())
        throw std::invalid_argument(
                "GridFunction::GridFunction(): "
                "functions from 'space' have a different number of "
                "components than 'function'");
    if (mode != APPROXIMATE && mode != INTERPOLATE)
        throw std::invalid_argument(
                "GridFunction::GridFunction(): "
                "'mode' must be either APPROXIMATE or INTERPOLATE");

    bool isBarycentricSpace = (space->isBarycentric() || dualSpace->isBarycentric());
    if (isBarycentricSpace) {
        m_space = space->barycentricSpace(space);
        m_dualSpace = dualSpace->barycentricSpace(dualSpace);
    }
    if (m_space->grid() != m_dualSpace->grid())
            throw std::invalid_argument(
                    "GridFunction::GridFunction(): "
                    "space and dualSpace must be defined on the same grid");

    if (mode == APPROXIMATE)
        setProjections(m_dualSpace,
                       *calculateProjections(*context, function, *m_dualSpace));
    else // mode == INTERPOLATE
        setCoefficients(interpolate(function, *m_space));
}

// Deprecated constructors

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
        const arma::Col<ResultType>& data,
        DataType dataType)
{
    bool isBarycentricSpace=false;
    if (space && space->isBarycentric())
        isBarycentricSpace = true;

    if (dualSpace && dualSpace->isBarycentric())
        isBarycentricSpace = true;

    shared_ptr<const Space<BasisFunctionType> > newSpace(space);
    shared_ptr<const Space<BasisFunctionType> > newDualSpace(dualSpace);
    if (isBarycentricSpace){
        newSpace = space->barycentricSpace(space);
        if (dualSpace)
            newDualSpace = dualSpace->barycentricSpace(dualSpace);
    }

    if (dataType == COEFFICIENTS) {
        if (newDualSpace && newSpace->grid() != newDualSpace->grid())
            throw std::invalid_argument(
                    "GridFunction::GridFunction(): "
                    "space and dualSpace must be defined on the same grid");
        initializeFromCoefficients(context, newSpace, data);
        m_dualSpace = newDualSpace;
    } else if (dataType == PROJECTIONS) {
        initializeFromProjections(context, newSpace, newDualSpace, data);
        m_dualSpace = newDualSpace;
    } else
        throw std::invalid_argument(
                "GridFunction::GridFunction(): invalid dataType");
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>::GridFunction(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
        const arma::Col<ResultType>& coefficients,
        const arma::Col<ResultType>& projections)
{
    bool isBarycentricSpace=false;
    if (space && space->isBarycentric())
        isBarycentricSpace = true;

    if (dualSpace && dualSpace->isBarycentric())
        isBarycentricSpace = true;

    shared_ptr<const Space<BasisFunctionType> > newSpace(space);
    shared_ptr<const Space<BasisFunctionType> > newDualSpace(dualSpace);
    if (isBarycentricSpace){
        newSpace = space->barycentricSpace(space);
        if (dualSpace)
            newDualSpace = dualSpace->barycentricSpace(dualSpace);
    }


    // We ignore the vector of projections.
    initializeFromCoefficients(context, newSpace, coefficients);
    if (newDualSpace && newSpace->grid() != newDualSpace->grid())
        throw std::invalid_argument(
                "GridFunction::GridFunction(): "
                "space and dualSpace must be defined on the same grid");
    m_dualSpace = newDualSpace;
}

// Member functions

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::initializeFromCoefficients(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const arma::Col<ResultType>& coefficients)
{
    if (!context)
        throw std::invalid_argument(
                "GridFunction::initializeFromCoefficients(): context must not be null");
    if (!space)
        throw std::invalid_argument(
                "GridFunction::initializeFromCoefficients(): space must not be null");
    if (coefficients.n_rows != space->globalDofCount())
        throw std::invalid_argument(
                "GridFunction::initializeFromCoefficients(): "
                "the coefficients vector has incorrect length");
    m_context = context;
    m_space = space;
    setCoefficients(coefficients);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::initializeFromProjections(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& space,
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace,
        const arma::Col<ResultType>& projections)
{
    if (!context)
        throw std::invalid_argument(
                "GridFunction::initializeFromProjections(): context must not be null");
    if (!space)
        throw std::invalid_argument(
                "GridFunction::initializeFromProjections(): space must not be null");
    if (!dualSpace)
        throw std::invalid_argument(
                "GridFunction::initializeFromProjections(): dualSpace must not be null");
    bool isBarycentricSpace = (space->isBarycentric() || dualSpace->isBarycentric());
    if (isBarycentricSpace) {
        m_space = space->barycentricSpace(space);
        m_dualSpace = dualSpace->barycentricSpace(dualSpace);
    }
    else {
        m_space = space;
        m_dualSpace = dualSpace;
    }
    if (m_space->grid() != m_dualSpace->grid())
            throw std::invalid_argument(
                    "GridFunction::initializeFromProjections(): "
                    "space and dualSpace must be defined on the same grid");

    if (projections.n_rows != dualSpace->globalDofCount())
        throw std::invalid_argument(
                "GridFunction::initializeFromProjections(): "
                "the projections vector has incorrect length");
    m_context = context;
    setProjections(m_dualSpace, projections);
}

template <typename BasisFunctionType, typename ResultType>
bool GridFunction<BasisFunctionType, ResultType>::isInitialized() const
{
    return m_space;
}

template <typename BasisFunctionType, typename ResultType>
bool GridFunction<BasisFunctionType, ResultType>::
wasInitializedFromCoefficients() const
{
    return m_wasInitializedFromCoefficients;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Grid> GridFunction<BasisFunctionType, ResultType>::grid() const
{
    if (!m_space)
        throw std::runtime_error("GridFunction::grid() must not be called "
                                 "on an uninitialized GridFunction object");

    return m_space->grid();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType> >
GridFunction<BasisFunctionType, ResultType>::space() const
{
    return m_space;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType> >
GridFunction<BasisFunctionType, ResultType>::dualSpace() const
{
    return m_dualSpace;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Context<BasisFunctionType, ResultType> >
GridFunction<BasisFunctionType, ResultType>::context() const
{
    return m_context;
}

template <typename BasisFunctionType, typename ResultType>
int GridFunction<BasisFunctionType, ResultType>::componentCount() const
{
    if (!m_space)
        throw std::runtime_error("GridFunction::componentCount() must not be called "
                                 "on an uninitialized GridFunction object");

    return m_space->codomainDimension();
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType,ResultType> GridFunction<BasisFunctionType,ResultType>::approximateInSpace(
        const shared_ptr<Space<BasisFunctionType> >& space ) const {

    BoundaryOperator<BasisFunctionType,ResultType> id =
            identityOperator<BasisFunctionType,ResultType> (m_context,m_space,space,space);
    return id*(*this);


}


template <typename BasisFunctionType, typename ResultType>
const arma::Col<ResultType>&
GridFunction<BasisFunctionType, ResultType>::coefficients() const
{
    if (!m_space || (!m_coefficients && !m_projections))
        throw std::runtime_error("GridFunction::coefficients() must not be called "
                                 "on an uninitialized GridFunction object");
    if (!m_coefficients)
        updateCoefficientsFromProjections();
    return *m_coefficients;
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setCoefficients(
        const arma::Col<ResultType>& coeffs)
{
    if (!m_space)
        throw std::runtime_error("GridFunction::setCoefficients() must not be called "
                                 "on an uninitialized GridFunction object");
    if (coeffs.n_rows != m_space->globalDofCount())
        throw std::invalid_argument(
                "GridFunction::setCoefficients(): dimension of the provided "
                "vector does not match the number of global DOFs in the primal space");
    m_coefficients.reset(new arma::Col<ResultType>(coeffs));
    m_projections.reset();
    m_wasInitializedFromCoefficients = true;
}

template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType>
GridFunction<BasisFunctionType, ResultType>::projections(
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace_) const
{
    if (!m_space)
        throw std::runtime_error("GridFunction::projections() must not be called "
                                 "on an uninitialized GridFunction object");
    if (m_space->grid() != dualSpace_->grid())
        if (!m_space->grid()->isBarycentricRepresentationOf(*dualSpace_->grid()) &&
                !dualSpace_->grid()->isBarycentricRepresentationOf(*m_space->grid()))
            throw std::invalid_argument(
                    "GridFunction::projections(): "
                    "space and dual space must be defined on the same grid");

    if (!m_coefficients && !m_projections)
        throw std::runtime_error("GridFunction::projections() must not be called "
                                 "on an uninitialized GridFunction object");
    if (!m_projections || !m_dualSpace->spaceIsCompatible(*dualSpace_))
        updateProjectionsFromCoefficients(dualSpace_);
    return *m_projections;
}

// Deprecated version
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType>
GridFunction<BasisFunctionType, ResultType>::projections(
        const Space<BasisFunctionType>& dualSpace_) const
{
    if (!m_dualSpace)
        throw std::runtime_error(
                "You must provide the dualSpace_ argument in the call to "
                "GridFunction::projections() if you did not specify the "
                "dual space when constructing the GridFunction.");
    arma::Col<ResultType> result = projections(
                make_shared_from_ref(dualSpace_));
    // update the coefficients because we can't guarantee
    // dualSpace_ will stay alive until they are needed
    coefficients();
    return result;
}

// Deprecated version
template <typename BasisFunctionType, typename ResultType>
arma::Col<ResultType>
GridFunction<BasisFunctionType, ResultType>::projections() const
{
    if (!m_dualSpace)
        throw std::runtime_error(
                "You must provide the dualSpace_ argument in the call to "
                "GridFunction::projections() if you did not specify the "
                "dual space when constructing the GridFunction.");
    return projections(m_dualSpace);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setProjections(
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace_,
        const arma::Col<ResultType>& projects)
{
    if (!m_space)
        throw std::runtime_error("GridFunction::setProjections() must not be called "
                                 "on an uninitialized GridFunction object");
    if (m_space->grid() != dualSpace_->grid())
        throw std::invalid_argument(
                "GridFunction::setProjections(): "
                "space and dual space must be defined on the same grid");
    if (projects.n_rows != dualSpace_->globalDofCount())
        throw std::invalid_argument(
                "GridFunction::setProjections(): dimension of the provided "
                "vector does not match the number of global DOFs in the dual space");

    m_projections = boost::make_shared<arma::Col<ResultType> >(projects);
    m_dualSpace = dualSpace_;
    m_coefficients.reset();
    m_wasInitializedFromCoefficients = false;
}

// Deprecated version
template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setProjections(
        const Space<BasisFunctionType>& dualSpace_,
        const arma::Col<ResultType>& projects)
{
    setProjections(make_shared_from_ref(dualSpace_), projects);
    // update the coefficients because we can't guarantee
    // dualSpace_ will stay alive until they are needed
    coefficients();
}

// Deprecated version
template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::setProjections(
        const arma::Col<ResultType>& projects)
{
    if (!m_dualSpace)
        throw std::runtime_error(
                "You must provide the dualSpace_ argument in the call to "
                "GridFunction::setProjections() if you did not specify the "
                "dual space when constructing the GridFunction.");
    setProjections(m_dualSpace, projects);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::
updateProjectionsFromCoefficients(
        const shared_ptr<const Space<BasisFunctionType> >& dualSpace_) const
{
    // This should have been checked beforehand, this is after all a private
    // function. So omit these checks in release mode.
    assert(isInitialized());
    assert(dualSpace_);
    assert(m_coefficients);

    // Calculate the mass matrix
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    BoundaryOp id = identityOperator(
                m_context, m_space, m_space, dualSpace_);

    shared_ptr<arma::Col<ResultType> > newProjections(
                new arma::Col<ResultType>(dualSpace_->globalDofCount()));
    id.weakForm()->apply(NO_TRANSPOSE, *m_coefficients, *newProjections,
                         static_cast<ResultType>(1.),
                         static_cast<ResultType>(0.));
    m_projections = newProjections;
    m_dualSpace = dualSpace_;
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::
updateCoefficientsFromProjections() const
{
    // This should have been checked beforehand, this is after all a private
    // function. So omit these checks in release mode.
    assert(isInitialized());
    assert(m_projections);
    assert(m_dualSpace);

    // Calculate the (pseudo)inverse mass matrix
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    BoundaryOp id = identityOperator(
                m_context, m_space, m_space, m_dualSpace);
    BoundaryOp pinvId = pseudoinverse(id);

    shared_ptr<arma::Col<ResultType> > newCoefficients(
                new arma::Col<ResultType>(m_space->globalDofCount()));
    pinvId.weakForm()->apply(
                NO_TRANSPOSE, *m_projections, *newCoefficients,
                static_cast<ResultType>(1.), static_cast<ResultType>(0.));
    m_coefficients = newCoefficients;
}

template <typename BasisFunctionType, typename ResultType>
typename GridFunction<BasisFunctionType, ResultType>::MagnitudeType
GridFunction<BasisFunctionType, ResultType>::L2Norm() const
{
    // The L^2 norm is given by
    //   sqrt(u^\dagger M u),
    // where u is the coefficient vector and M the mass matrix of m_space

    if (!m_space)
        throw std::runtime_error("GridFunction::L2_Norm() must not be called "
                                 "on an uninitialized GridFunction object");
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;

    // Get the vector of coefficients
    const arma::Col<ResultType>& coeffs = coefficients();

    // Calculate the mass matrix
    BoundaryOp id = identityOperator(m_context, m_space, m_space, m_space);
    shared_ptr<const DiscreteBoundaryOperator<ResultType> > massMatrix =
            id.weakForm();

    arma::Col<ResultType> product(coeffs.n_rows);
    massMatrix->apply(NO_TRANSPOSE, coeffs, product, 1., 0.);
    ResultType result = arma::cdot(coeffs, product);
    if (fabs(imagPart(result)) > 1000. * std::numeric_limits<MagnitudeType>::epsilon())
        std::cout << "Warning: squared L2Norm has non-negligible imaginary part: "
                  << imagPart(result) << std::endl;
    return sqrt(realPart(result));
}

// Redundant, in fact -- can be obtained directly from Space
template <typename BasisFunctionType, typename ResultType>
const Fiber::Shapeset<BasisFunctionType>&
GridFunction<BasisFunctionType, ResultType>::shapeset(
        const Entity<0>& element) const
{
    BOOST_ASSERT_MSG(m_space, "GridFunction::shapeset() must not be "
                     "called on an uninitialized GridFunction object");
    return m_space->shapeset(element);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::getLocalCoefficients(
        const Entity<0>& element, std::vector<ResultType>& coeffs) const
{
    BOOST_ASSERT_MSG(m_space, "GridFunction::getLocalCoefficients() must not be "
                     "called on an uninitialized GridFunction object");
    std::vector<GlobalDofIndex> gdofIndices;
    std::vector<BasisFunctionType> ldofWeights;
    m_space->getGlobalDofs(element, gdofIndices, ldofWeights);
    const int gdofCount = gdofIndices.size();
    coeffs.resize(gdofCount);
    const arma::Col<ResultType>& globalCoefficients = coefficients();
    for (int i = 0; i < gdofCount; ++i) {
        int gdof = gdofIndices[i];
        if (gdof >= 0)
            coeffs[i] = globalCoefficients(gdof) * ldofWeights[i];
        else
            coeffs[i] = 0.;
    }
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::exportToVtk(
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType) const
{
    Bempp::exportToVtk(*this, dataType, dataLabel,
                       fileNamesBase, filesPath, outputType);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::evaluateAtSpecialPoints(
        VtkWriter::DataType dataType, arma::Mat<ResultType>& result_) const
{
    arma::Mat<CoordinateType> points;
    evaluateAtSpecialPoints(dataType, points, result_);
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::evaluateAtSpecialPoints(
        VtkWriter::DataType dataType, arma::Mat<CoordinateType>& points,
        arma::Mat<ResultType>& values) const
{
    if (!m_space)
        throw std::runtime_error("GridFunction::evaluateAtSpecialPoints() must "
                                 "not be called on an uninitialized GridFunction object");
    if (dataType != VtkWriter::CELL_DATA && dataType != VtkWriter::VERTEX_DATA)
        throw std::invalid_argument("GridFunction::evaluateAtSpecialPoints(): "
                                    "invalid data type");

    const GridView& view = m_space->gridView();
    const int gridDim = m_space->gridDimension();
    const int worldDim = m_space->worldDimension();
    const int elementCodim = 0;
    const int vertexCodim = gridDim;
    const int nComponents = componentCount();


    const size_t elementCount = view.entityCount(elementCodim);
    const size_t vertexCount = view.entityCount(vertexCodim);

    values.set_size(nComponents,
                     dataType == VtkWriter::CELL_DATA ? elementCount : vertexCount);
    values.fill(0.);
    points.set_size(worldDim, values.n_cols);

    // Number of elements contributing to each column in result
    // (this will be greater than 1 for VERTEX_DATA)
    std::vector<int> multiplicities(vertexCount);
    std::fill(multiplicities.begin(), multiplicities.end(), 0);

    // Gather geometric data
    Fiber::RawGridGeometry<CoordinateType> rawGeometry(
                gridDim, m_space->worldDimension());
    view.getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData(), rawGeometry.domainIndices());

    // Make geometry factory
    shared_ptr<const Grid> grid = m_space->grid();
    std::auto_ptr<GeometryFactory> geometryFactory =
            grid->elementGeometryFactory();
    std::auto_ptr<typename GeometryFactory::Geometry> geometry(
                geometryFactory->make());
    Fiber::GeometricalData<CoordinateType> geomData;

    // For each element, get its shapeset and corner count (this is sufficient
    // to identify its geometry) as well as its local coefficients
    typedef std::pair<const Fiber::Shapeset<BasisFunctionType>*, int>
            ShapesetAndCornerCount;
    typedef std::vector<ShapesetAndCornerCount> ShapesetAndCornerCountVector;
    ShapesetAndCornerCountVector basesAndCornerCounts(elementCount);
    std::vector<std::vector<ResultType> > localCoefficients(elementCount);
    {
        const Mapper& mapper = view.elementMapper();
        std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
        for (size_t e = 0; e < elementCount; ++e) {
            const Entity<0>& element = it->entity();
            const int elementIndex = mapper.entityIndex(element);
            basesAndCornerCounts[elementIndex] = ShapesetAndCornerCount(
                        &m_space->shapeset(element),
                        rawGeometry.elementCornerCount(elementIndex));
            getLocalCoefficients(element, localCoefficients[elementIndex]);
            it->next();
        }
    }

    typedef std::set<ShapesetAndCornerCount> ShapesetAndCornerCountSet;
    ShapesetAndCornerCountSet uniqueShapesetsAndCornerCounts(
                basesAndCornerCounts.begin(), basesAndCornerCounts.end());

    // Find out which basis data need to be calculated
    size_t basisDeps = 0, geomDeps = Fiber::GLOBALS;
    // Find out which geometrical data need to be calculated, in addition
    // to those needed by the kernel
    const Fiber::CollectionOfShapesetTransformations<CoordinateType>& transformations =
            m_space->basisFunctionValue();
    assert(nComponents == transformations.resultDimension(0));
    transformations.addDependencies(basisDeps, geomDeps);

    // Loop over unique combinations of basis and element corner count
    typedef typename ShapesetAndCornerCountSet::const_iterator
            BasisAndCornerCountSetConstIt;
    for (BasisAndCornerCountSetConstIt it = uniqueShapesetsAndCornerCounts.begin();
         it != uniqueShapesetsAndCornerCounts.end(); ++it) {
        const ShapesetAndCornerCount& activeBasisAndCornerCount = *it;
        const Fiber::Shapeset<BasisFunctionType>& activeShapeset =
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
        activeShapeset.evaluate(basisDeps, local, ALL_DOFS, basisData);

        Fiber::BasisData<ResultType> functionData;
        if (basisDeps & Fiber::VALUES)
            functionData.values.set_size(basisData.values.extent(0),
                                         1, // just one function
                                         basisData.values.extent(2));
        if (basisDeps & Fiber::DERIVATIVES)
            functionData.derivatives.set_size(basisData.derivatives.extent(0),
                                              basisData.derivatives.extent(1),
                                              1, // just one function
                                              basisData.derivatives.extent(3));
        Fiber::CollectionOf3dArrays<ResultType> functionValues;

        // Loop over elements and process those that use the active shapeset
        for (size_t e = 0; e < elementCount; ++e) {
            if (basesAndCornerCounts[e].first != &activeShapeset)
                continue;

            // Local coefficients of the argument in the current element
            const std::vector<ResultType>& activeLocalCoefficients =
                    localCoefficients[e];

            // Calculate the function's values and/or derivatives
            // at the requested points in the current element
            if (basisDeps & Fiber::VALUES) {
                std::fill(functionData.values.begin(),
                          functionData.values.end(), 0.);
                for (size_t point = 0; point < basisData.values.extent(2); ++point)
                    for (size_t dim = 0; dim < basisData.values.extent(0); ++dim)
                        for (size_t fun = 0; fun < basisData.values.extent(1); ++fun)
                            functionData.values(dim, 0, point) +=
                                    basisData.values(dim, fun, point) *
                                    activeLocalCoefficients[fun];
            }
            if (basisDeps & Fiber::DERIVATIVES) {
                std::fill(functionData.derivatives.begin(),
                          functionData.derivatives.end(), 0.);
                for (size_t point = 0; point < basisData.derivatives.extent(3); ++point)
                    for (size_t dim = 0; dim < basisData.derivatives.extent(1); ++dim)
                        for (size_t comp = 0; comp < basisData.derivatives.extent(0); ++comp)
                            for (size_t fun = 0; fun < basisData.derivatives.extent(2); ++fun)
                                functionData.derivatives(comp, dim, 0, point) +=
                                        basisData.derivatives(comp, dim, fun, point) *
                                        activeLocalCoefficients[fun];
            }

            // Get geometrical data
            rawGeometry.setupGeometry(e, *geometry);
            geometry->getData(geomDeps, local, geomData);
            if (geomDeps & Fiber::DOMAIN_INDEX)
                geomData.domainIndex = rawGeometry.domainIndex(e);

            transformations.evaluate(functionData, geomData, functionValues);
            assert(functionValues[0].extent(1) == 1); // one function

            if (dataType == VtkWriter::CELL_DATA) {
                for (int dim = 0; dim < nComponents; ++dim)
                    values(dim, e) = functionValues[0](     // array index
                                                       dim, // component
                                                       0,   // function index
                                                       0);  // point index
                for (int dim = 0; dim < worldDim; ++dim)
                    points(dim, e) = geomData.globals(dim, 0);
            }
            else { // VERTEX_DATA
                // Add the calculated values to the columns of the result array
                // corresponding to the active element's vertices
                for (int c = 0; c < activeCornerCount; ++c) {
                    int vertexIndex = rawGeometry.elementCornerIndices()(c, e);
                    for (int dim = 0; dim < nComponents; ++dim)
                        values(dim, vertexIndex) += functionValues[0](dim, 0, c);
                    ++multiplicities[vertexIndex];
                }
                for (int c = 0; c < activeCornerCount; ++c) {
                    int vertexIndex = rawGeometry.elementCornerIndices()(c, e);
                    for (int dim = 0; dim < worldDim; ++dim)
                        points(dim, vertexIndex) = geomData.globals(dim, c);
                }
            }
        } // end of loop over elements
    } // end of loop over unique combinations of shapeset and corner count

    // Take average of the vertex values obtained in each of the adjacent elements
    if (dataType == VtkWriter::VERTEX_DATA)
        for (size_t v = 0; v < vertexCount; ++v)
            values.col(v) /= multiplicities[v];
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::evaluate(
        const Entity<0>& element,
        const arma::Mat<CoordinateType>& local,
        arma::Mat<ResultType>& values) const
{
    if (local.n_rows != m_space->grid()->dim())
        throw std::invalid_argument("evaluate(): points in 'local' have an "
                                    "invalid number of coordinates");

    const int nComponents = componentCount();
    // Find out which basis data need to be calculated
    size_t basisDeps = 0, geomDeps = 0;
    // Find out which geometrical data need to be calculated,
    const Fiber::CollectionOfShapesetTransformations<CoordinateType>& transformations =
        m_space->basisFunctionValue();
    assert(transformations.transformationCount() == 1);
    assert(nComponents == transformations.resultDimension(0));
    transformations.addDependencies(basisDeps, geomDeps);

    // Get basis data
    const Fiber::Shapeset<BasisFunctionType>& shapeset = m_space->shapeset(element);
    Fiber::BasisData<BasisFunctionType> basisData;
    shapeset.evaluate(basisDeps, local, ALL_DOFS, basisData);
    // Get geometrical data
    Fiber::GeometricalData<CoordinateType> geomData;
    element.geometry().getData(geomDeps, local, geomData);
    values.set_size(nComponents, local.n_cols);
    // Get shape function values
    Fiber::CollectionOf3dArrays<BasisFunctionType> functionValues;
    transformations.evaluate(basisData, geomData, functionValues);
    assert(functionValues.size() == 1);

    // Get local coefficients
    std::vector<ResultType> localCoefficients(shapeset.size());
    getLocalCoefficients(element, localCoefficients);
    assert(localCoefficients.size() == shapeset.size());

    // Calculate grid function values
    values.set_size(functionValues[0].extent(0), local.n_cols);
    values.fill(static_cast<ResultType>(0.));
    for (size_t p = 0; p < functionValues[0].extent(2); ++p)
        for (size_t f = 0; f < functionValues[0].extent(1); ++f)
            for (size_t dim = 0; dim < functionValues[0].extent(0); ++dim)
                values(dim, p) += functionValues[0](dim, f, p) * localCoefficients[f];
}

template <typename BasisFunctionType, typename ResultType>
void GridFunction<BasisFunctionType, ResultType>::exportToGmsh(
        const char* dataLabel, const char* fileName) const
{
    Bempp::exportToGmsh(*this, dataLabel, fileName);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g)
{
    return g;
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g)
{
    return static_cast<ResultType>(-1.) * g;
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


template <typename BasisFunctionType, typename ResultType>
void exportToVtk(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        VtkWriter::DataType dataType,
        const char* dataLabel,
        const char* fileNamesBase, const char* filesPath,
        VtkWriter::OutputType outputType)
{
    shared_ptr<const Space<BasisFunctionType> > space = gridFunction.space();
    if (!space)
        throw std::runtime_error("exportToVtk(): gridFunction must not be "
                                 "an uninitialized GridFunction object");
    arma::Mat<ResultType> data;
    gridFunction.evaluateAtSpecialPoints(dataType, data);

    std::auto_ptr<GridView> view = space->grid()->leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();

    exportSingleDataSetToVtk(*vtkWriter, data, dataType, dataLabel,
                             fileNamesBase, filesPath, outputType);
}


BEMPP_GCC_DIAG_OFF(deprecated-declarations);

// Redundant, in fact -- can be obtained directly from Space
template <typename BasisFunctionType, typename ResultType>
const Fiber::Basis<BasisFunctionType>&
GridFunction<BasisFunctionType, ResultType>::basis(
        const Entity<0>& element) const
{
    BOOST_ASSERT_MSG(m_space, "GridFunction::basis() must not be "
                     "called on an uninitialized GridFunction object");
    return m_space->basis(element);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator+(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2)
{
    if (g1.space() != g2.space())
        throw std::runtime_error("GridFunction::operator+(): spaces don't match");
    if (g1.wasInitializedFromCoefficients() ||
            g2.wasInitializedFromCoefficients() ||
            g1.dualSpace() != g2.dualSpace())
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(), // arbitrary choice...
                    g1.space(),
                    g1.coefficients() + g2.coefficients());
    else {
        shared_ptr<const Space<BasisFunctionType> > dualSpace = g1.dualSpace();
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(), // arbitrary choice...
                    g1.space(),
                    dualSpace,
                    g1.projections(dualSpace) + g2.projections(dualSpace));
    }
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator-(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const GridFunction<BasisFunctionType, ResultType>& g2)
{
    if (g1.space() != g2.space())
        throw std::runtime_error("GridFunction::operator-(): spaces don't match");
    // For the sake of old-style code (with dualSpace stored in the GridFunction),
    // we try to provide a sensible dual space to the composite GridFunction.
    if (g1.wasInitializedFromCoefficients() ||
            g2.wasInitializedFromCoefficients() ||
            g1.dualSpace() != g2.dualSpace())
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(), // arbitrary choice...
                    g1.space(),
                    g1.coefficients() - g2.coefficients());
    else {
        shared_ptr<const Space<BasisFunctionType> > dualSpace = g1.dualSpace();
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(), // arbitrary choice...
                    g1.space(),
                    dualSpace,
                    g1.projections(dualSpace) - g2.projections(dualSpace));
    }
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const GridFunction<BasisFunctionType, ResultType>& g1,
        const ScalarType& scalar)
{
    if (g1.wasInitializedFromCoefficients())
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(),
                    g1.space(),
                    static_cast<ResultType>(scalar) * g1.coefficients());
    else {
        shared_ptr<const Space<BasisFunctionType> > dualSpace = g1.dualSpace();
        return GridFunction<BasisFunctionType, ResultType>(
                    g1.context(),
                    g1.space(),
                    dualSpace,
                    static_cast<ResultType>(scalar) * g1.projections(dualSpace));
    }
}

BEMPP_GCC_DIAG_ON(deprecated-declarations);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template GridFunction<BASIS, RESULT> operator+( \
    const GridFunction<BASIS, RESULT>& op); \
    template GridFunction<BASIS, RESULT> operator-( \
    const GridFunction<BASIS, RESULT>& op); \
    template GridFunction<BASIS, RESULT> operator+( \
    const GridFunction<BASIS, RESULT>& op1, \
    const GridFunction<BASIS, RESULT>& op2); \
    template GridFunction<BASIS, RESULT> operator-( \
    const GridFunction<BASIS, RESULT>& op1, \
    const GridFunction<BASIS, RESULT>& op2); \
    template void exportToVtk( \
    const GridFunction<BASIS, RESULT>& gridFunction, \
    VtkWriter::DataType dataType, \
    const char* dataLabel, \
    const char* fileNamesBase, const char* filesPath, \
    VtkWriter::OutputType outputType)
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
// INSTANTIATE_FREE_FUNCTIONS_FOR_REAL_BASIS(float);
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
// INSTANTIATE_FREE_FUNCTIONS_FOR_COMPLEX_BASIS(std::complex<float>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(
        double, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, double);
// INSTANTIATE_FREE_FUNCTIONS_FOR_REAL_BASIS(double);
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
// INSTANTIATE_FREE_FUNCTIONS_FOR_COMPLEX_BASIS(std::complex<double>);
#endif

} // namespace Bempp
