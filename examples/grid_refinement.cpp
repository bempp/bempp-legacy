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
#include "grid/reverse_element_mapper.hpp"
#include "grid/grid_view.hpp"
#include "grid/geometry.hpp"
#include "grid/entity.hpp"
#include "grid/entity_pointer.hpp"
#include "grid/mapper.hpp"



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

     using namespace Bempp;

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

      const char* meshFile = "/Users/betcke/local/bempp/development-debug/bempp/examples/meshes/sphere-h-0.4.msh";
      GridParameters params;
      params.topology = GridParameters::TRIANGULAR;
      shared_ptr<Grid> grid = GridFactory::importGmshGrid(params, meshFile);
      shared_ptr<Grid> baryGrid = grid->barycentricGrid();

      PiecewiseLinearContinuousScalarSpace<BFT> plins(grid);
      PiecewiseConstantDualMeshScalarSpaceBarycentric<BFT> plinsd(grid);

      std::vector<GlobalDofIndex> globalDofIndex;
      globalDofIndex.push_back(0);

      std::vector<std::vector<LocalDof> > localDofs;

      plinsd.global2localDofs(globalDofIndex,localDofs);

      const std::vector<LocalDof>& zeroDof = localDofs[0];
      std::auto_ptr<GridView> view = baryGrid->leafView();
      std::auto_ptr<GridView> viewCoarse = baryGrid->levelView(0);
      const ReverseElementMapper& reverseMapper = view->reverseElementMapper();
      const Mapper& mapper = viewCoarse->elementMapper();




      for (int i=0;i<zeroDof.size();++i){
          std::cout << "("<< zeroDof[i].entityIndex<<","<<zeroDof[i].dofIndex<<")"<< std::endl;
          const EntityPointer<0>& elementPointer = reverseMapper.entityPointer(zeroDof[i].entityIndex);
          const Entity<0>& element = elementPointer.entity();
          const Geometry& geom = element.geometry();
          arma::Mat<double> corners;
          geom.getCorners(corners);
          std::cout << corners;
          std::auto_ptr<EntityPointer<0> > father = element.father();
          const Entity<0>& fatherElement = father->entity();
          std::cout << "Father element: " << mapper.entityIndex(fatherElement) << std::endl;
      }

      std::cout << "Linear space" << std::endl;

      plins.global2localDofs(globalDofIndex,localDofs);
      const std::vector<LocalDof>& zeroDof2 = localDofs[0];
      for (int i=0;i<zeroDof2.size();++i){
          std::cout << "("<< zeroDof2[i].entityIndex<<","<<zeroDof2[i].dofIndex<<")"<<std::endl;
      }

}
