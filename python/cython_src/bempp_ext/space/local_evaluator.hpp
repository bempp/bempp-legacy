#ifndef BEMPP_EXT_LOCAL_EVALUATOR_HPP
#define BEMPP_EXT_LOCAL_EVALUATOR_HPP

#include "bempp/common/common.hpp"
#include "bempp/common/eigen_support.hpp"
#include "bempp/common/types.hpp"
#include "bempp/space/space.hpp"
#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/basis.hpp"
#include "bempp/fiber/basis_data.hpp"
#include "bempp/fiber/collection_of_basis_transformations.hpp"
#include "bempp/grid/geometry.hpp"
#include "bempp/grid/grid.hpp"
#include "bempp/grid/entity.hpp"

namespace Bempp {


template <typename ValueType>
Matrix<ValueType> evaluateLocalBasis(const Space<double>& space, const Entity<0>& element, 
        const Matrix<double>& local, const Vector<ValueType>& localCoefficients)
{


  if (local.rows() != space.grid()->dim())
    throw std::invalid_argument("evaluate(): points in 'local' have an "
                                "invalid number of coordinates");

  const int nComponents = space.codomainDimension();
  Matrix<ValueType> values(nComponents, local.cols());
  values.setZero();

  // Find out which basis data need to be calculated
  size_t basisDeps = 0, geomDeps = 0;
  // Find out which geometrical data need to be calculated,
  const Fiber::CollectionOfShapesetTransformations<double>
      &transformations = space.basisFunctionValue();
  transformations.addDependencies(basisDeps, geomDeps);

  // Get basis data
  const Fiber::Shapeset<double> &shapeset =
      space.shapeset(element);
  Fiber::BasisData<double> basisData;
  shapeset.evaluate(basisDeps, local, ALL_DOFS, basisData);
  // Get geometrical data
  Fiber::GeometricalData<double> geomData;
  element.geometry().getData(geomDeps, local, geomData);
  // Get shape function values
  Fiber::CollectionOf3dArrays<double> functionValues;
  transformations.evaluate(basisData, geomData, functionValues);

  // Calculate grid function values
  for (size_t p = 0; p < functionValues[0].extent(2); ++p)
    for (size_t f = 0; f < functionValues[0].extent(1); ++f)
      for (size_t dim = 0; dim < functionValues[0].extent(0); ++dim)
        values(dim, p) += functionValues[0](dim, f, p) * localCoefficients(f);

  return values;
}

}

#endif
