// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_integration_manager_hpp
#define fiber_integration_manager_hpp

#include "array_2d.hpp"
#include "types.hpp"
#include <vector>
#include <stdexcept>

namespace Fiber
{

template <typename ValueType> class Basis;
template <typename ValueType, typename GeometryImp> class DoubleIntegrator;

// Some stream of consciousness...

// Additional quadrature order specific to the Helmholtz equation (dependence
// on k*h) will need to be handled via a kernel-specific code (there will be a
// particular quadrature selector for Helmholtz kernels, probably derived from
// the standard Sauter-Schwab quadrature selector).

// The quadrature manager needs to know the polynomial order of an element
// and the increment in the integration order resulting from mesh curvature.

/** \brief Integration manager interface.

  This class template is used as a base class for all integration manager
  implementations. It uses the Barton-Nackman trick to ensure conformity to the
  interface.

  An integration manager provides methods that choose appropriate quadrature
  rules (weights and points) for integrals occurring in boundary-element
  matrices. Factors influencing the choice of quadrature rule include the
  regularity properties of a kernel, maximum polynomial order of basis
  functions on an element, element volume and (for double integrals)
  element-element distance.
 */
template <typename ValueType, typename GeometryImp>
class IntegrationManager
{
public:
//    /** \brief Number of dimensions over which integration takes place. */
//    int elementDimension() const {
//        asImp().elementDimension();
//    }

    void getTestKernelTrialIntegrators(
            CallVariant callVariant,
            const std::vector<const GeometryImp*>& geometriesA,
            const GeometryImp& geometryB,
            const std::vector<const Basis<ValueType>*>& basesA,
            const Basis<ValueType>& basisB,
            std::vector<const DoubleIntegrator<ValueType, GeometryImp>*>& integrators) {
        const int elementCount = geometriesA.size();
        if (basesA.size() != elementCount)
            throw std::invalid_argument("IntegrationManager::"
                                        "getTestKernelTrialIntegrators(): "
                                        "incompatible argument lengths");
        integrators.resize(elementCount);
        if (callVariant = TEST_TRIAL)
            for (int i = 0; i < elementCount; ++i)
                integrators[i] = &testKernelTrialIntegrator(
                            *geometriesA[i], geometryB, *basesA[i], basisB);
        else
            for (int i = 0; i < elementCount; ++i)
                integrators[i] = &testKernelTrialIntegrator(
                            geometryB, *geometriesA[i], basisB, *basesA[i]);
    }

    void getTestKernelTrialIntegrators(
            const std::vector<const GeometryImp*>& testGeometries,
            const std::vector<const GeometryImp*>& trialGeometries,
            const std::vector<const Basis<ValueType>*>& testBases,
            const std::vector<const Basis<ValueType>*>& trialBases,
            Array2D<const DoubleIntegrator<ValueType, GeometryImp>*>& integrators) {
        const int testElementCount = testGeometries.size();
        const int trialElementCount = trialGeometries.size();
        if (testBases.size() != testElementCount ||
                trialBases.size() != trialElementCount)
            throw std::invalid_argument("IntegrationManager::"
                                        "getTestKernelTrialIntegrators(): "
                                        "incompatible argument lengths");
        integrators.set_size(testElementCount, trialElementCount);
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                integrators(testIndex, trialIndex) = &testKernelTrialIntegrator(
                            *testGeometries[testIndex],
                            *trialGeometries[trialIndex],
                            *testBases[testIndex],
                            *trialBases[trialIndex]);
    }

    // non-const because it might create a new integrator internally
    virtual const DoubleIntegrator<ValueType, GeometryImp>& testKernelTrialIntegrator(
            const GeometryImp& testGeometry,
            const GeometryImp& trialGeometry,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis) = 0;
};

} // namespace Fiber

#endif
