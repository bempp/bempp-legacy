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

#ifndef fiber_double_integrator_hpp
#define fiber_double_integrator_hpp

namespace Fiber
{

/** \brief Integration over pairs of elements. */
class DoubleIntegrator
{
public:
    // Note: trialFamily or testFamily can very well contain only one function
    // if we are interested in one local DOF only.

    /** \brief Integrate test functions times kernel times trial functions on a
      set of test and trial elements.

      Below, elementDimension is the value returned by the method dimension()
      of the geometries (the dimension of all elements must be the same).

      \param[in]  pointCount
                  Number of quadrature points.

      \param[in]  localTestQuadPoints

                  Pointer to a Fortran-ordered 2D array of dimensions
                  (elementDimension, pointCount) storing the local coordinates
                  of quadrature points on the test elements.

      \param[in]  localTrialQuadPoints

                  Pointer to a Fortran-ordered 2D array of dimensions
                  (elementDimension, pointCount) storing the local coordinates
                  of quadrature points on the test elements.

      \param[in]  quadWeigths

                  Pointer to an array of length pointCount storing the weights
                  of quadrature points.

      \param[in]  testGeometryCount

                  Number of test elements.

      \param[in]  testGeometries

                  Pointer to a testGeometryCount-long array of constant
                  pointers to objects implementing the Geometry interface and
                  representing the mappings from local to global coordinates
                  on test elements.

      \param[in]  trialGeometryCount

                  Number of trial elements.

      \param[in]  trialGeometries

                  Pointer to a trialGeometryCount-long array of constant
                  pointers to objects implementing the Geometry interface and
                  representing the mappings from local to global coordinates
                  on trial elements.

      \param[in]  testFamily

                  An object implementing the FunctionFamily interface and
                  representing a collection of test functions defined on the
                  test elements.

      \param[in]  trialFamily

                  An object implementing the FunctionFamily interface and
                  representing a collection of trial functions defined on the
                  trial element.

      \param[in]  kernel

                  An object implementing the Kernel interface and
                  representing the kernel of an integral operator.

      \param[out] result

                  Pointer to a preallocated Fortran-ordered 3D array of
                  dimensions (testFamily.size(), testGeometries.size(),
                  trialFamily.size()), which will receive the calculated
                  integrals. */
    template <typename ValueType,
              typename CoordinateType,
              typename GeometryImp,
              typename TestFunctionFamilyImp,
              typename TrialFunctionFamilyImp,
              typename KernelImp>
    static void integrate(
            int pointCount,
            const CoordinateType* localTestQuadPoints,
            const CoordinateType* localTrialQuadPoints,
            const ValueType* quadWeights,
            int testGeometryCount,
            const GeometryImp** testGeometries,
            int trialGeometryCount,
            const GeometryImp** trialGeometries,
            TestFunctionFamilyImp& testFamily,
            TrialFunctionFamilyImp& trialFamily,
            const KernelImp& kernel,
            ValueType* result);
};

} // namespace Fiber

#include "double_integrator_imp.hpp"

#endif
