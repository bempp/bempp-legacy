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
#include <cmath>
#include <iostream>
#include <memory>

#include "meshes.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_linear_operator.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "assembly/linear_operator_superposition.hpp"

using namespace Bempp;
using std::cout;
using std::endl;

int main(void){


    const MeshVariant meshVariant = CUBE_12; // SPHERE_2590;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);

    //PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    PiecewiseConstantScalarSpace<double> space(*grid);

    space.assignDofs();

    AssemblyOptions assemblyOptions;
    AcaOptions acaOptions;
    assemblyOptions.switchToAca(acaOptions);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    typedef std::auto_ptr<ElementaryLinearOperator<double> > LinearOperatorPtr;
    IdentityOperator<double> ident(space,space);
    SingleLayerPotential3D<double> slp(space,space);
    LinearOperatorSuperposition<double> idSlp=2.*ident+slp;


    std::auto_ptr<DiscreteLinearOperator<double> > discreteId =
            ident.assembleWeakForm(factory, assemblyOptions);
    std::auto_ptr<DiscreteLinearOperator<double> > discreteSlp =
            slp.assembleWeakForm(factory,assemblyOptions);
    std::auto_ptr<DiscreteLinearOperator<double> > discreteIdSlp =
            idSlp.assembleWeakForm(factory,assemblyOptions);


    std::cout << "Identity" << std::endl;
    discreteId->dump();

    std::cout << "SLP" << std::endl;
    discreteSlp->dump();

    std::cout << "Sum" << std::endl;
    discreteIdSlp->dump();

    return 0;
}
