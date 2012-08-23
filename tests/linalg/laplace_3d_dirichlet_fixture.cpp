#include "laplace_3d_dirichlet_fixture.hpp"

#include "assembly/context.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"
#include "fiber/explicit_instantiation.hpp"
#include "grid/grid_factory.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

namespace Bempp
{

template <typename BFT, typename RT>
Laplace3dDirichletFixture<BFT, RT>::Laplace3dDirichletFixture(
    SpaceType dirichletDataDomain,
    SpaceType neumannDataDomain,
    SpaceType range,
    SpaceType dualToRange)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    grid = GridFactory::importGmshGrid(
        params, "../examples/meshes/cube-12-reoriented.msh");

    shared_ptr<Space<BFT> > pwiseConstants(
        new PiecewiseConstantScalarSpace<BFT>(*grid));
    shared_ptr<Space<BFT> > pwiseLinears(
        new PiecewiseLinearContinuousScalarSpace<BFT>(*grid));
    pwiseConstants->assignDofs();
    pwiseLinears->assignDofs();

    AssemblyOptions assemblyOptions;
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy( 
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));
        
    BoundaryOperator<BFT, RT> slpOp = laplace3dSingleLayerBoundaryOperator<BFT, RT>(
        context, 
        neumannDataDomain == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants, 
        range == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants,
        dualToRange == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants);
    BoundaryOperator<BFT, RT> dlpOp = laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
        context, 
        dirichletDataDomain == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants, 
        range == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants,
        dualToRange == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants);
    BoundaryOperator<BFT, RT> id = identityOperator<BFT, RT>(
        context, 
        dirichletDataDomain == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants, 
        range == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants,
        dualToRange == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants);
        
    lhsOp = slpOp;
    BoundaryOperator<BFT, RT> rhsOp = -0.5 * id + dlpOp;
        
    GridFunction<BFT, RT> u(
        context,
        dirichletDataDomain == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants, 
        dirichletDataDomain == PIECEWISE_LINEARS ? pwiseLinears : pwiseConstants, 
        surfaceNormalIndependentFunction(UnitFunctor<RT>()));
    rhs = rhsOp * u;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dDirichletFixture);

} // namespace Bempp
