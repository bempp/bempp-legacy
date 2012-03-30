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
#include <memory>

#include "meshes.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "assembly/assembly_options.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"
#include "assembly/double_layer_potential_3d.hpp"
#include "assembly/adjoint_double_layer_potential_3d.hpp"
#include "assembly/hypersingular_operator_3d.hpp"
#include "assembly/linear_operator_superposition.hpp"
#include "assembly/elementary_linear_operator.hpp"
#include "assembly/discrete_scalar_valued_linear_operator_superposition.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/geometry.hpp"

#include "fiber/accuracy_options.hpp"
#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

using namespace Bempp;

/**
 * Represents the (Laplace) Calderon Projection on a given mesh.
 */
class Calderon {
    std::auto_ptr<DiscreteScalarValuedLinearOperator<double> > A, B, C, D;
public:
	Calderon(Grid* grid){
		PiecewiseLinearContinuousScalarSpace<double> space1(*grid);
		PiecewiseLinearContinuousScalarSpace<double> space2(*grid);
		space1.assignDofs();
		space2.assignDofs();

		AssemblyOptions assemblyOptions;
		assemblyOptions.switchToDense();

		Fiber::AccuracyOptions accuracyOptions; // default
		Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
				factory(accuracyOptions);

        typedef std::auto_ptr<ElementaryLinearOperator<double> > ELOP;

        SingleLayerPotential3D<double> slp;
        ELOP dblneg(new DoubleLayerPotential3D<double>());
        dblneg->scale(-1.0);
        ELOP adj(new AdjointDoubleLayerPotential3D<double>());
        HypersingularOperator3D<double> hyp;

        ELOP halfid1(new IdentityOperator<double>());
        halfid1->scale(0.5);
        ELOP halfid2(new IdentityOperator<double>());
        halfid2->scale(0.5);

        typedef LinearOperatorSuperposition<double> LOSd;

        A = LOSd(halfid1, dblneg).assembleWeakForm(space1, space1, factory, assemblyOptions);
        B = slp.assembleWeakForm(space2, space1, factory, assemblyOptions);
        C = hyp.assembleWeakForm(space1, space2, factory, assemblyOptions);
        D = LOSd(halfid2, adj).assembleWeakForm(space2, space2, factory, assemblyOptions);
	}


    arma::Mat<double> getMatrix(){
    	return arma::join_cols(arma::join_rows(A->asMatrix(), B->asMatrix()),
    						   arma::join_rows(C->asMatrix(), D->asMatrix()));
    }


};

int main(){
    const MeshVariant meshVariant = CUBE_384;
    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    Calderon c(grid.get());
    arma::Mat<double> m = c.getMatrix();
    arma::Mat<double> m2 = m*m;

    std::cout<<"M"<<std::endl<<c.getMatrix().submat(0,0,20,20);

    std::cout<<"M*M"<<std::endl<<m2.submat(0,0,20,20);

    //    std::cout<<"M*M*M"<<std::endl<<(m*m*m).submat(0,0,20,20);

}




