// Copyright (C) 2011-2012 by the Bem++ Authors
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

#include "../common/common.hpp"

// Keep IDEs happy
#include "default_local_assembler_for_local_operators_on_surfaces.hpp"

#include "numerical_test_trial_integrator.hpp"
#include "quadrature_descriptor_selector_for_local_operators.hpp"
#include "single_quadrature_rule_family.hpp"

#include <boost/tuple/tuple_comparison.hpp>

namespace Fiber {

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
DefaultLocalAssemblerForLocalOperatorsOnSurfaces<BasisFunctionType, ResultType,
                                                 GeometryFactory>::
    DefaultLocalAssemblerForLocalOperatorsOnSurfaces(
        const shared_ptr<const GeometryFactory> &geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<
            CoordinateType>> &testTransformations,
        const shared_ptr<const CollectionOfShapesetTransformations<
            CoordinateType>> &trialTransformations,
        const shared_ptr<const TestTrialIntegral<BasisFunctionType, ResultType>>
            &integral,
        const shared_ptr<const OpenClHandler> &openClHandler,
        const shared_ptr<const QuadratureDescriptorSelectorForLocalOperators<
            CoordinateType>> &quadDescSelector,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
            &quadRuleFamily)
    : m_geometryFactory(geometryFactory), m_rawGeometry(rawGeometry),
      m_testShapesets(testShapesets), m_trialShapesets(trialShapesets),
      m_testTransformations(testTransformations),
      m_trialTransformations(trialTransformations), m_integral(integral),
      m_openClHandler(openClHandler), m_quadDescSelector(quadDescSelector),
      m_quadRuleFamily(quadRuleFamily) {}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForLocalOperatorsOnSurfaces<
    BasisFunctionType, ResultType, GeometryFactory>::
    evaluateLocalWeakForms(const std::vector<int> &elementIndices,
                           std::vector<Matrix<ResultType>> &result) {
  typedef TestTrialIntegrator<BasisFunctionType, ResultType> Integrator;
  typedef Shapeset<BasisFunctionType> Shapeset;

  const int elementCount = elementIndices.size();
  result.resize(elementCount);

  // Select integrator for each element
  typedef boost::tuples::tuple<const Integrator *, const Shapeset *,
                               const Shapeset *> QuadVariant;
  std::vector<QuadVariant> quadVariants(elementCount);
  for (int i = 0; i < elementCount; ++i) {
    const Integrator *integrator = &selectIntegrator(elementIndices[i]);
    quadVariants[i] =
        QuadVariant(integrator, (*m_testShapesets)[elementIndices[i]],
                    (*m_trialShapesets)[elementIndices[i]]);
  }

  // Integration will proceed in batches of test elements having the same
  // "quadrature variant", i.e. integrator and shapesets

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<int> activeElementIndices;
  activeElementIndices.reserve(elementCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *it->template get<0>();
    const Shapeset &activeTestShapeset = *it->template get<1>();
    const Shapeset &activeTrialShapeset = *it->template get<2>();

    // Find all the test elements for which quadrature should proceed
    // according to the current quadrature variant
    activeElementIndices.clear();
    for (int e = 0; e < elementCount; ++e)
      if (quadVariants[e] == activeQuadVariant)
        activeElementIndices.push_back(elementIndices[e]);

    // Integrate!
    std::vector<Matrix<ResultType>> localResult;
    activeIntegrator.integrate(activeElementIndices, activeTestShapeset,
                               activeTrialShapeset, localResult);

    // Distribute the just calculated integrals into the result array
    // that will be returned to caller
    int i = 0;
    for (int e = 0; e < elementCount; ++e)
      if (quadVariants[e] == activeQuadVariant)
        result[e] = localResult[i++];
  }
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
const TestTrialIntegrator<BasisFunctionType, ResultType> &
DefaultLocalAssemblerForLocalOperatorsOnSurfaces<
    BasisFunctionType, ResultType,
    GeometryFactory>::selectIntegrator(int elementIndex) {
  SingleQuadratureDescriptor desc =
      m_quadDescSelector->quadratureDescriptor(elementIndex);
  return getIntegrator(desc);
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
const TestTrialIntegrator<BasisFunctionType, ResultType> &
DefaultLocalAssemblerForLocalOperatorsOnSurfaces<
    BasisFunctionType, ResultType,
    GeometryFactory>::getIntegrator(const SingleQuadratureDescriptor &desc) {
  typename IntegratorMap::iterator it = m_testTrialIntegrators.find(desc);
  if (it != m_testTrialIntegrators.end()) {
    //            std::cout << "getIntegrator(: " << index << "): integrator
    // found" << std::endl;
    return *it->second;
  }
  //        std::cout << "getIntegrator(: " << index << "): integrator not
  // found" << std::endl;

  // Integrator doesn't exist yet and must be created.
  Matrix<CoordinateType> points;
  std::vector<CoordinateType> weights;
  m_quadRuleFamily->fillQuadraturePointsAndWeights(desc, points, weights);

  typedef NumericalTestTrialIntegrator<BasisFunctionType, ResultType,
                                       GeometryFactory> Integrator;
  std::unique_ptr<TestTrialIntegrator<BasisFunctionType, ResultType>>
      integrator(new Integrator(points, weights, *m_geometryFactory,
                                *m_rawGeometry, *m_testTransformations,
                                *m_trialTransformations, *m_integral,
                                *m_openClHandler));

  SingleQuadratureDescriptor key(desc);
  return *m_testTrialIntegrators.insert(key, integrator.release())
              .first->second;
}

} // namespace Fiber
