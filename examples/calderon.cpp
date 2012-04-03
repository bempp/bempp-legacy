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
using namespace arma;

/**
 * Represents the (Laplace) Calderon Projection on a given mesh.
 */
class Calderon {
    std::auto_ptr<DiscreteScalarValuedLinearOperator<double> > slp, adj, hyp, dbl, id12, id21;
public:
	Calderon(Grid* grid){
		PiecewiseLinearContinuousScalarSpace<double> space1(*grid);
//		PiecewiseLinearContinuousScalarSpace<double> space2(*grid);
		PiecewiseConstantScalarSpace<double> space2(*grid);
		space1.assignDofs();
		space2.assignDofs();

		AssemblyOptions assemblyOptions;
		assemblyOptions.switchToDense();

		Fiber::AccuracyOptions accuracyOptions; // default
		accuracyOptions.regular.orderIncrement = 0;
		accuracyOptions.singular.orderIncrement = 0;
		Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
				factory(accuracyOptions);

        dbl = DoubleLayerPotential3D<double>().assembleWeakForm(space1, space2, factory, assemblyOptions);
		slp = SingleLayerPotential3D<double>().assembleWeakForm(space2, space2, factory, assemblyOptions);
		hyp = HypersingularOperator3D<double>().assembleWeakForm(space1, space1, factory, assemblyOptions);
		adj = AdjointDoubleLayerPotential3D<double>().assembleWeakForm(space2, space1, factory, assemblyOptions);

		id21 = IdentityOperator<double>().assembleWeakForm(space2,space1,factory, assemblyOptions);
        id12 = IdentityOperator<double>().assembleWeakForm(space1,space2,factory, assemblyOptions);
	}


    Mat<double> getMatrix(){
//    	cout<<id12->asMatrix()<<id21->asMatrix();
    	mat ID12inv = pinv(id12->asMatrix());
    	mat ID21inv = pinv(id21->asMatrix());
    	mat DBL = ID21inv * dbl->asMatrix().t();
    	mat SLP = ID21inv * slp->asMatrix();
    	mat HYP = ID12inv * hyp->asMatrix();
    	mat ADJ = ID12inv * adj->asMatrix().t();
    	mat I1 = eye(SLP.n_rows, SLP.n_rows);
    	mat I2 = eye(HYP.n_rows, HYP.n_rows);
		mat M = join_cols(join_rows(I1 * 0.5 - DBL, SLP),
						  join_rows(HYP, I2* 0.5 + ADJ));
		return M;
    }

};

int main(){
    const MeshVariant meshVariant = SPHERE_614;
    std::auto_ptr<Grid> grid = loadMesh(meshVariant);
    Calderon c(grid.get());
    Mat<double> m = c.getMatrix();
    Mat<double> m2 = m*m;

    m.save("matrix.csv", csv_ascii);

//    std::cout<<"M"<<std::endl<<m;
//    std::cout<<"M*M"<<std::endl<<m2;

    mat U, V;
    vec S;
    std::cout<<svd(U,S,V,m - m2)<<std::endl;
    std::cout<<"Singular values of difference"<<std::endl<<S.t();
    std::cout<<svd(U,S,V,m)<<std::endl;
    std::cout<<"Singular values of M"<<std::endl<<S.t();


    //    std::cout<<"M"<<std::endl<<c.getMatrix().submat(0,0,20,20);
//
//    std::cout<<"M*M"<<std::endl<<m2.submat(0,0,20,20);

    //    std::cout<<"M*M*M"<<std::endl<<(m*m*m).submat(0,0,20,20);

}




