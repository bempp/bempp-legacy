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

#ifndef fiber_default_quadrature_descriptor_selector_for_potential_operators_hpp
#define fiber_default_quadrature_descriptor_selector_for_potential_operators_hpp

#include "quadrature_descriptor_selector_for_potential_operators.hpp"

#include "../common/shared_ptr.hpp"
#include "accuracy_options.hpp"
#include "scalar_traits.hpp"
#include "types.hpp"

namespace Fiber {

template <typename BasisFunctionType> class Shapeset;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType>
class DefaultLocalAssemblerForOperatorsOnSurfacesUtilities;

/** \brief Default implementation of a quadrature descriptor selector
 *  used during the evaluation of potentials.
 *
 *  The choice of quadrature rule accuracy can be influenced by the
 *  \p accuracyOptions parameter taken by the constructor. */
template <typename BasisFunctionType>
class DefaultQuadratureDescriptorSelectorForPotentialOperators
    : public QuadratureDescriptorSelectorForPotentialOperators<
          BasisFunctionType> {
  typedef QuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>
      Base;

public:
  typedef typename Base::CoordinateType CoordinateType;

  DefaultQuadratureDescriptorSelectorForPotentialOperators(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const AccuracyOptionsEx &accuracyOptions);

  virtual SingleQuadratureDescriptor
  quadratureDescriptor(const Vector<CoordinateType> &point,
                       int trialElementIndex,
                       CoordinateType nominalDistance) const;

  virtual SingleQuadratureDescriptor
  farFieldQuadratureDescriptor(const Shapeset<BasisFunctionType> &trialShapeset,
                               int trialElementCornerCount) const;

private:
  /** \cond PRIVATE */
  void precalculateElementSizesAndCenters();

  int order(const Vector<CoordinateType> &point, int trialElementIndex,
            CoordinateType nominalDistance) const;
  CoordinateType
  pointElementDistanceSquared(const Vector<CoordinateType> &point,
                              int trialElementIndex) const;

  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

  shared_ptr<const RawGridGeometry<CoordinateType>> m_rawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
      m_trialShapesets;
  AccuracyOptionsEx m_accuracyOptions;

  std::vector<CoordinateType> m_elementSizesSquared;
  Matrix<CoordinateType> m_elementCenters;
  CoordinateType m_averageElementSize;
  /** \endcond */
};

} // namespace Fiber

#endif
