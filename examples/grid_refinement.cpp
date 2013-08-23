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

#include "assembly/assembly_options.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"
#include "assembly/discrete_boundary_operator.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "linalg/default_iterative_solver.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_constant_dual_mesh_scalar_space_barycentric.hpp"
#include "space/piecewise_linear_continuous_scalar_space_barycentric.hpp"
#include "space/piecewise_linear_discontinuous_scalar_space_barycentric.hpp"

#include <iostream>
#include <fstream>

typedef double BFT; // basis function type
typedef double RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

class DirichletData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = 2 * x * z / (r * r * r * r * r) - y / (r * r * r);
    }
};

class ExactNeumannData
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
    typedef CT CoordinateType;

    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    // Number of components of the function's value
    int resultDimension() const { return 1; }

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        CoordinateType r = sqrt(point(0) * point(0) +
                point(1) * point(1) +
                point(2) * point(2));
        result(0) = -6 * x * z / (r * r * r * r * r * r) + 2 * y / (r * r * r * r);
    }
};

int main()
{
//    // Import symbols from namespace Bempp to the global namespace

//    using namespace Bempp;

//    // Load mesh

//    arma::Col<double> lowerLeft(2);
//    arma::Col<double> upperRight(2);
//    arma::Col<unsigned int> elemNumbers(2);

//    arma::Mat<double> corners(3,3);
//    arma::Mat<int> elements(3,1);

//    corners(0,0) = 0;
//    corners(1,0) = 0;
//    corners(2,0) = 0;

//    corners(0,1) = 1;
//    corners(1,1) = 0;
//    corners(2,1) = 0;

//    corners(0,2) = 0;
//    corners(1,2) = 1;
//    corners(2,2) = 1;

//    elements(0,0) = 0;
//    elements(1,0) = 1;
//    elements(2,0) = 2;


////    lowerLeft(0) = 0;
////    lowerLeft(1) = 0;
////    upperRight(0) = 1;
////    upperRight(1) = 1;
////    elemNumbers(0) = 1;
////    elemNumbers(1) = 1;

//    //const char* meshFile = "/Users/betcke/local/bempp/development-debug/bempp/examples/meshes/sphere-h-0.4.msh";
//    GridParameters params;
//    params.topology = GridParameters::TRIANGULAR;
//    //shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);
//    //shared_ptr<Grid> grid = GridFactory::createStructuredGrid(params,lowerLeft,upperRight,elemNumbers);
//    shared_ptr<Grid> grid = GridFactory::createGridFromConnectivityArrays(params,corners,elements);
//    //shared_ptr<Grid> grid2 = grid->barycentricGrid();

//    // Initialize the spaces

//    //PiecewiseConstantDualMeshScalarSpaceBarycentric<BFT> pconsts(grid2);
//    //PiecewiseConstantDualMeshScalarSpaceBarycentric<BFT> pconsts(grid);
//    //std::cout << "Dual Space end" << std::endl;

//    PiecewiseLinearContinuousScalarSpaceBarycentric<BFT> plins(grid);
//    PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BFT>plinsd(grid);




//    // Define the quadrature strategy

//    AccuracyOptions accuracyOptions;
//    // Increase by 2 the order of quadrature rule used to approximate
//    // integrals of regular functions on pairs on elements
//    accuracyOptions.doubleRegular.setRelativeQuadratureOrder(8);
//    accuracyOptions.doubleSingular.setRelativeQuadratureOrder(8);
//    // Increase by 2 the order of quadrature rule used to approximate
//    // integrals of regular functions on single elements
//    accuracyOptions.singleRegular.setRelativeQuadratureOrder(6);
//    NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);

//    // Specify the assembly method. We want to use ACA

//    AssemblyOptions assemblyOptions;
//    AcaOptions acaOptions; // Default parameters for ACA
//    assemblyOptions.switchToAcaMode(acaOptions);

//    // Create the assembly context

//    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

//    // Construct elementary operators

//    BoundaryOperator<BFT, RT> hypersing =
//            laplace3dHypersingularBoundaryOperator<BFT, RT>(
//                make_shared_from_ref(context),
//                make_shared_from_ref(plins),
//                make_shared_from_ref(plins),
//                make_shared_from_ref(plins));

////    BoundaryOperator<BFT, RT> slpOpd =
////            laplace3dSingleLayerBoundaryOperator<BFT, RT>(
////                make_shared_from_ref(context),
////                make_shared_from_ref(plinsd),
////                make_shared_from_ref(plinsd),
////                make_shared_from_ref(plinsd));


//    shared_ptr<const DiscreteBoundaryOperator<RT> > weak = hypersing.weakForm();
//    std::cout << "("<<weak->rowCount()<<","<<weak->columnCount()<<")"<<std::endl;
////    shared_ptr<const DiscreteBoundaryOperator<RT> > slpWeakd = slpOpd.weakForm();

////    arma::Mat<double> m = slpWeak->asMatrix();
////    arma::Mat<double> md = slpWeakd->asMatrix();

////    std::cout << m << " " << std::endl;
////    std::cout << md << " " << std::endl;


////    std::cout << "("<<slpWeak->rowCount()<<","<<slpWeak->columnCount()<<")"<<std::endl;
////    std::cout << "("<<slpWeakd->rowCount()<<","<<slpWeakd->columnCount()<<")"<<std::endl;
////    return 0;

    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasisBarycentric<BFT> Basis;
    Basis basis(Basis::TYPE1);
    arma::Mat<Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<BFT> data;
    basis.evaluate(Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    Fiber::_4dArray<BFT> expected(1, // component count
                                        2, //
                                        vertexCount,
                                        pointCount);
    std::fill(expected.begin(), expected.end(), 0.);

    expected(0, 0, 0, 0) = -2./3.;
    expected(0, 1, 0, 0) = -1./2.;
    expected(0, 0, 0, 1) = -2./3.;
    expected(0, 1, 0, 1) = -1./2.;
    expected(0, 0, 0, 2) = -2./3.;
    expected(0, 1, 0, 2) = -1./2.;

    expected(0, 0, 1, 0) = 1./3.;
    expected(0, 1, 1, 0) = 0.;
    expected(0, 0, 1, 1) = 1./3.;
    expected(0, 1, 1, 1) = 0.;
    expected(0, 0, 1, 2) = 1./3.;
    expected(0, 1, 1, 2) = 0.;

    expected(0, 0, 2, 0) = 1./3.;
    expected(0, 1, 2, 0) = 1./2;
    expected(0, 0, 2, 1) = 1./3.;
    expected(0, 1, 2, 1) = 1./2;
    expected(0, 0, 2, 2) = 1./3.;
    expected(0, 1, 2, 2) = 1./2;

    for (int i=0;i<expected.extent(0);++i)
        for (int j=0;j<expected.extent(1);++j)
            for (int k=0;k<expected.extent(2);++k)
                for (int l=0;l<expected.extent(3);++l){
                    std::cout << "("<<i<<","<<j<<","<<k<<","<<l<<"): "<< expected(i,j,k,l) << ", "<< data.derivatives(i,j,k,l)<< std::endl;
                }
}
