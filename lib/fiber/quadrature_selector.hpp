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

#ifndef fiber_quadrature_selector_hpp
#define fiber_quadrature_selector_hpp

namespace Fiber
{

///** \brief Basic parameters of a quadrature rule. */
//struct QuadratureRule
//{
//    /** \brief Unique identifier of the quadrature rule. */
//    int id;
//    /** \brief Number of quadrature points in this rule. */
//    int pointCount;
//    bool operator==(const DoubleQuadratureRule& other) const {
//        return (id == other.id &&
//                pointCount == other.pointCount);
//    }
//};

// Some stream of consciousness...

// Additional quadrature order specific to the Helmholtz equation (dependence
// on k*h) will need to be handled via a kernel-specific code (there will be a
// particular quadrature selector for Helmholtz kernels, probably derived from
// the standard Sauter-Schwab quadrature selector).
//
// What if an equation contains several very different kernels in an equation?
// Maybe each integral operator should provide a factory method generating *its
// own* well-adapted quadrature selector, with the user able to override this
// choice. All this doesn't seem to influence the basic code design and can be
// fine-tuned later.

typedef int QuadratureRule;

/** \brief Quadrature selector interface.

  This class template is used as a base class for all quadrature selector
  implementations. It uses the Barton-Nackman trick to ensure conformity to the
  interface.

  A quadrature selector provides methods that choose appropriate quadrature
  rules (weights and points) for integrals occurring in boundary-element
  matrices. Factors influencing the choice of quadrature rule include the
  regularity properties of a kernel, maximum polynomial order of basis
  functions on an element, element volume and (for double integrals)
  element-element distance.

  \tparam ValueType      Type used to represent quadrature weights.
  \tparam CoordinateType Type used to represent components of quadrature points.
  \tparam KernelImp      Type of an implementation of the quadrature selector
                         interface. */
template <typename ValueType, typename CoordinateType,
          typename QuadratureSelectorImp>
class QuadratureSelector
{
public:
    /** \brief Number of dimensions over which integration takes place. */
    int elementDimension() const {
        asImp().elementDimension();
    }

    /** \brief Select appropriate quadratures rule for integration over pairs
      of (test and trial) elements.

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
      \param[out] quadRules
                  Pointer to a preallocated Fortran-ordered 2D array of
                  dimensions (testGeometryCount, trialGeometryCount), in which
                  will be stored the identifiers of the quadrature rules
                  selected for each pair of elements. */
    template <typename GeometryImp>
    void selectDoubleQuadratureRules(
            int testGeometryCount,
            const GeometryImp** testGeometries,
            int trialGeometryCount,
            const GeometryImp** trialGeometries,
            QuadratureRule* quadRules) const {
        asImp().selectDoubleQuadratureRules(testGeometryCount, testGeometries,
                                            trialGeometryCount, trialGeometries,
                                            quadRules);
    }

    /** \brief Get number of points in the given double quadrature rule. */
    int doubleQuadratureRulePointCount(QuadratureRule quadRule) const {
        return asImp().doubleQuadratureRulePointCount(quadRule);
    }

    /** \brief Get weights and points corresponding to the given double
      quadrature.

      Below, pointCount is the value returned by
      doubleQuadratureRulePointCount(quadRule).

      \param[in]  quadRule
                  Quadrature rule id.
      \param[out] testQuadPoints
                  Pointer to a preallocated Fortran-ordered 2D array of
                  dimensions (elementDimension(), pointCount),
                  in which will be stored the coordinates of quadrature points
                  on the test elements.
      \param[out] trialQuadPoints
                  Pointer to a preallocated Fortran-ordered 2D array of
                  dimensions (elementDimension(), pointCount),
                  in which will be stored the coordinates of quadrature points
                  on the trial elements.
      \param[out] quadWeights
                  Pointer to a preallocated array of length pointCount,
                  in which will be stored the quadrature weights. */
    void doubleQuadratureRulePointsAndWeights(
            QuadratureRule quadRule,
            CoordinateType* testQuadPoints,
            CoordinateType* trialQuadPoints,
            ValueType* quadWeights) const {
        asImp().doubleQuadratureRulePointsAndWeights(
                    quadRule, testQuadPoints, trialQuadPoints, quadWeights);
    }

private:
    /**  Barton-Nackman trick */
    QuadratureSelectorImp& asImp() {
        return static_cast<QuadratureSelectorImp&>(*this);
    }

    /**  Barton-Nackman trick */
    const QuadratureSelectorImp& asImp() const {
        return static_cast<const QuadratureSelectorImp&>(*this);
    }
};

} // namespace Fiber

#endif
