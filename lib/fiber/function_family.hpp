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

#ifndef fiber_function_family_hpp
#define fiber_function_family_hpp

namespace Fiber
{

/** \brief Function family interface.

  This class template is used as a base class for all implementations of
  function families. It uses the Barton-Nackman trick to ensure conformity to
  the interface.

  An implementation of this interface represents a set of one or more possibly
  vector-valued shape functions defined on a reference element. Evaluation of
  these functions proceeds in two steps: first, the coordinates of the
  evaluation points must be passed to setEvaluationPoints(), and then the
  function values may be retrieved by calling evaluate(). The reason for this
  two-step procedure is that the functions will usually depend on the basis
  functions (defined on the reference element), which can be precalculated and
  cached.

  \tparam ValueType      Type used to represent components of the function values
                         (e.g. double or std::complex<float>).
  \tparam CoordinateType Type used to represent components of the evaluation
                         point coordinates (e.g. float or double).
  \tparam FamilyImp      Type of an implementation of the function family
                         interface. */
template <typename ValueType, typename CoordinateType, typename FamilyImp>
class FunctionFamily
{
public:
    /** \brief Dimension of the reference element on which the functions are defined. */
    int domainDimension() const {
        return asImp().domainDimension();
    }

    /** \brief Dimension of the space in which the element is embedded. */
    int worldDimension() const {
        return asImp().worldDimension();
    }

    /** \brief Number of components of the function values. */
    int codomainDimension() const {
        return asImp().codomainDimension();
    }

    /** \brief Number of functions belonging to the family. */
    int size() const {
        return asImp().size();
    }

    /** \brief True if the functions depend on the (pseudo)inverse of the
        transposed Jacobian matrix of the coordinate mapping from the reference
        to the physical element. */
    bool needsJacobianInverseTransposed() const {
        return asImp().needsJacobianInverseTransposed();
    }

    /** \brief Set points at which the functions will be evaluated.

      \param[in] pointCount Number of points.
      \param[in] pointCount Pointer to a Fortran-ordered 2D array of dimensions
                            (domainDimension(), pointCount) storing local
                            coordinates of points at which the functions will be
                            evaluated. */
    void setEvaluationPoints(int pointCount,
                             const CoordinateType* evaluationPoints) {
        asImp().setEvaluationPoints(pointCount, evaluationPoints);
    }

    /** \brief Number of points supplied in the previous call to
      setEvaluationPoints() or zero if that function has not been called yet. */
    int evaluationPointCount() const {
        return asImp().getEvaluationPoints();
    }

    /** \brief Evaluate the functions at the points previously supplied to
      setEvaluationPoints().

      \param[in]  jacobianInversesTransposed      
                  Pointer to a Fortran-ordered 3D array of dimensions
                  (worldDimension(), domainDimension(), evaluationPointCount())
                  storing the (pseudo)inverses of the transposed Jacobian
                  matrices of the coordinate mapping from the reference to the
                  physical element at each of the evaluation points.
                  transformation described by this geometry at the given
                  points. May be NULL if needsJacobianInverseTransposed()
                  returns false.
      \param[out] result
                  Pointer to a preallocated Fortran-ordered 3D array of
                  dimensions (codomainDimension(), size(),
                  evaluationPointCount()), in which will be stored the
                  components of the functions belonging to the family,
                  evaluated at the points given in the previous call to
                  setEvaluationPoints(). */
    void evaluate(const CoordinateType* jacobianInversesTransposed,
                  ValueType* result) const {
        asImp().evaluate(jacobianInversesTransposed, result);
    }

private:
    /**  Barton-Nackman trick */
    FamilyImp& asImp() {
        return static_cast<FamilyImp&>(*this);
    }

    /**  Barton-Nackman trick */
    const FamilyImp& asImp() const {
        return static_cast<const FamilyImp&>(*this);
    }
};

} // namespace Fiber

#endif
