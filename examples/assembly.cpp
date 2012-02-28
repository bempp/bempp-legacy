// Copyright (C) 2011 by the BEM++ Authors
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

#include <armadillo>
#include <iostream>
#include <memory> // auto_ptr

#include "assembly/assembly_options.hpp"
#include "assembly/single_layer_potential_3d.hpp"

#include "assembly/discrete_scalar_valued_linear_operator.hpp"
#include "fiber/standard_integration_manager_factory_2d.hpp"

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"
#include "grid/index_set.hpp"
#include "grid/vtk_writer.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"

using namespace Bempp;
using std::cout;
using std::endl;

/**
    A script for rudimentary testing of the single-layer-potential operator.
  */
int main()
{
    const char MESH_FNAME[] = "simple_mesh_9_elements.msh";

    // Import the grid
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    std::auto_ptr<Grid> grid(GridFactory::importGmshGrid(params, std::string(MESH_FNAME),
                             true, // verbose
                             false)); // insertBoundarySegments

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    assemblyOptions.mode = ASSEMBLY_MODE_DENSE;

    Fiber::OpenClOptions openClOptions;
    openClOptions.useOpenCl = false;
    Fiber::StandardIntegrationManagerFactory2D<double, GeometryFactory> factory(openClOptions);

    SingleLayerPotential3D<double> op;
    std::auto_ptr<DiscreteScalarValuedLinearOperator<double> > result =
            op.assembleWeakForm(space, space, factory, assemblyOptions);
    result->dump();
}



