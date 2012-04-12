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
#include "assembly/linear_operator_superposition.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/geometry.hpp"
#include "grid/geometry_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"

using namespace Bempp;

inline bool approxEqual(double x, double y)
{
    return fabs(x - y) / fabs((x + y) / 2.) < 1e-3;
}

inline bool approxZero(double x)
{
    return fabs(x) < 1e-6;
}

/**
    Demonstration of the usage of composite (superposition) operators
  */
int main()
{
    const MeshVariant meshVariant = CUBE_12_REORIENTED;

    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    dumpElementList(grid.get());

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    //assemblyOptions.switchToDense();

    AcaOptions acaOptions;
    acaOptions.eps = 1e-4;
    acaOptions.maximumRank = 10000;
    acaOptions.minimumBlockSize = 2;
    acaOptions.eta = 0.8;
    acaOptions.recompress = true;
    assemblyOptions.switchToAca(acaOptions);

    // assemblyOptions.switchToSparse();

    assemblyOptions.switchToTbb(1);
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::YES);

    // Fiber::OpenClOptions openClOptions;
    // assemblyOptions.switchToOpenCl(openClOptions);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    typedef std::auto_ptr<ElementaryLinearOperator<double> > LinearOperatorPtr;
    LinearOperatorPtr halfId(new IdentityOperator<double>);
    halfId->scale(0.5);
    LinearOperatorPtr slp(new SingleLayerPotential3D<double>);


    assemblyOptions.switchToSparse();
    std::auto_ptr<DiscreteLinearOperator<double> > resultHalfId =
            halfId->assembleWeakForm(space, space, factory, assemblyOptions);
    assemblyOptions.switchToDense();
    std::auto_ptr<DiscreteLinearOperator<double> > resultSlp =
            slp->assembleWeakForm(space, space, factory, assemblyOptions);

    assemblyOptions.switchToAca(acaOptions);

    LinearOperatorSuperposition<double> superposition(halfId, slp);
    std::auto_ptr<DiscreteLinearOperator<double> > resultSuperposition =
            superposition.assembleWeakForm(space, space, factory, assemblyOptions);

    std::cout << "\nHalfId:\n"
              << resultHalfId->asMatrix()
              << std::endl;

    std::cout << "\nSlp:\n"
              << resultSlp->asMatrix()
              << std::endl;

    std::cout << "\nSum of matrices generated separately:\n"
              << resultHalfId->asMatrix() + resultSlp->asMatrix()
              << std::endl;

    std::cout << "\nMatrix generated in a single run:\n"
              << resultSuperposition->asMatrix()
              << std::endl;
}



