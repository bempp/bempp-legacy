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

#ifndef fiber_default_quadrature_descriptor_selector_for_local_operators_hpp
#define fiber_default_quadrature_descriptor_selector_for_local_operators_hpp

#include "quadrature_descriptor_selector_for_local_operators.hpp"

#include "../common/shared_ptr.hpp"
#include "accuracy_options.hpp"
#include "scalar_traits.hpp"

namespace Fiber {

template <typename BasisFunctionType> class Shapeset;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType>
class DefaultLocalAssemblerForOperatorsOnSurfacesUtilities;

/** \brief Default implementation of a quadrature descriptor selector
 *  used during the discretization of local boundary operators. */
template <typename BasisFunctionType>
class DefaultQuadratureDescriptorSelectorForLocalOperators
    : public QuadratureDescriptorSelectorForLocalOperators<
          typename ScalarTraits<BasisFunctionType>::RealType> {
public:
  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  /** \note: \p accuracyOptions is currently unused. */
  DefaultQuadratureDescriptorSelectorForLocalOperators(
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const AccuracyOptionsEx &accuracyOptions);

  virtual SingleQuadratureDescriptor
  quadratureDescriptor(int elementIndex) const;

private:
  /** \cond PRIVATE */
  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

  shared_ptr<const RawGridGeometry<CoordinateType>> m_rawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_testShapesets;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_trialShapesets;
  AccuracyOptionsEx m_accuracyOptions;
  /** \endcond */
};

} // namespace Fiber

#endif
