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

#ifndef fiber_quadrature_descriptor_selector_factory_hpp
#define fiber_quadrature_descriptor_selector_factory_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "scalar_traits.hpp"

#include <vector>

namespace Fiber {

template <typename BasisFunctionType> class Shapeset;
template <typename CoordinateType> class RawGridGeometry;
template <typename CoordinateType>
class QuadratureDescriptorSelectorForGridFunctions;
template <typename CoordinateType>
class QuadratureDescriptorSelectorForIntegralOperators;
template <typename CoordinateType>
class QuadratureDescriptorSelectorForLocalOperators;
template <typename CoordinateType>
class QuadratureDescriptorSelectorForPotentialOperators;

/** \ingroup quadrature
 *  \brief Builder of quadrature descriptor selectors.
 *
 *  This class serves as a factory of quadrature descriptor selectors
 *  for particular grids. Quadrature descriptor selectors are objects
 *  that determine the type and degree of accuracy of quadrature rules
 *  needed to evaluate particular integrals.
 *
 *  The parameters of the \c make...() member functions provide
 *  information about the grids and shape functions used during the
 *  discretization of operators and functions. These data may be
 *  stored in the quadrature descriptor selectors being constructed
 *  and subsequently used to determine the parameters of quadrature
 *  rules. */
template <typename BasisFunctionType>
class QuadratureDescriptorSelectorFactory {
public:
  /** \brief Type used to represent coordinates. */
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  /** \brief Destructor. */
  virtual ~QuadratureDescriptorSelectorFactory() {}

  /** \brief Create a quadrature descriptor selector used during the
   *  discretization of functions.
   *
   *  The test functions used during the discretization live on the
   *  grid described by the \p rawGeometry object. The \p
   *  testShapesets vector contains pointers to the sets of shape
   *  functions defined on all the elements of this grid. */
  virtual shared_ptr<
      QuadratureDescriptorSelectorForGridFunctions<CoordinateType>>
  makeQuadratureDescriptorSelectorForGridFunctions(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &testShapesets) const = 0;

  /** \brief Create a quadrature descriptor selector used during
   *  the discretization of the weak form of boundary integral operators.
   *
   *  The test and trial functions used during the discretization
   *  live on the grids described by the \p testRawGeometry and \p
   *  trialRawGeometry objects, respectively. The \p testShapesets
   *  and \p trialShapesets vectors contain pointers to the sets of
   *  shape functions defined on all the elements of these grids. */
  virtual shared_ptr<
      QuadratureDescriptorSelectorForIntegralOperators<CoordinateType>>
  makeQuadratureDescriptorSelectorForIntegralOperators(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &trialShapesets) const = 0;

  /** \brief Create a quadrature descriptor selector used during
   *  the discretization of the weak form of local boundary operators.
   *
   *  The test and trial functions used during the discretization
   *  live on the grid described by the \p rawGeometry object. The
   *  \p testShapesets and \p trialShapesets vectors contain
   *  pointers to the sets of test and trial shape functions defined
   *  on all the elements of this grid. */
  virtual shared_ptr<
      QuadratureDescriptorSelectorForLocalOperators<CoordinateType>>
  makeQuadratureDescriptorSelectorForLocalOperators(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &trialShapesets) const = 0;

  /** \brief Create a quadrature descriptor selector used during
   *  the evaluation of potentials.
   *
   *  The trial functions used during the evaluation live on the
   *  grid described by the \p rawGeometry object. The \p
   *  trialShapesets vector contains pointers to the sets of shape
   *  functions defined on all the elements of this grid. */
  virtual shared_ptr<
      QuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>>
  makeQuadratureDescriptorSelectorForPotentialOperators(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
          &trialShapesets) const = 0;
};

} // namespace Fiber

#endif
