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

#ifndef fiber_kernel_hpp
#define fiber_kernel_hpp

// #include "barton_nackman.hpp"

namespace Fiber
{

/** \brief Kernel interface.

  This class template is used as a base class for all kernel implementations.
  It uses the Barton-Nackman trick to ensure conformity to the interface.

  A kernel implementation represents a particular integral operator's kernel,
  i.e. a tensor-valued function \f$K(x, y)\f$, where \f$x\f$ and \f$y\f$ are
  called test and trial points, respectively.

  \tparam ValueType      Type used to represent components of the kernel tensor
                         (e.g. double or std::complex<float>).
  \tparam CoordinateType Type used to represent components of the test and trial
                         point coordinates (e.g. float or double).
  \tparam KernelImp      Type of an implementation of the kernel interface. */
template <typename ValueType, typename CoordinateType,
          typename KernelImp>
class Kernel
{
public:
    /** \brief Number of components of kernel's arguments \f$x\f$ and \f$y\f$. */
    int worldDimension() const {
        return asImp().worldDimension();
    }

    /** \brief Number of rows of the tensor \f$K\f$. */
    int codomainDimension() const {
        return asImp().codomainDimension();
    }

    /** \brief Number of columns of the tensor \f$K\f$. */
    int domainDimension() const {
        return asImp().domainDimension();
    }

    /** \brief True if the kernel depends on the vector normal to the surface at
        test point. */
    bool needsTestNormal() const {
        return asImp().needsTestNormal();
    }

    /** \brief True if the kernel depends on the vector normal to the surface at
        trial point. */
    bool needsTrialNormal() const {
        return asImp().needsTrialNormal();
    }

    /** \brief Evaluate kernel at prescribed pairs of test and trial points.

      \param[in]  pointCount
                  Number of points
      \param[in]  testPoints
                  Pointer to a Fortran-ordered 2D array of dimensions
                  (worldDimension(), pointCount) storing the test point
                  coordinates.
      \param[in]  trialPoints
                  Pointer to a Fortran-ordered 2D array of dimensions
                  (worldDimension(), pointCount) storing the trial point
                  coordinates.
      \param[in]  testNormals
                  Pointer to a Fortran-ordered 2D array of dimensions
                  (worldDimension(), pointCount) storing the components of
                  external unit vectors normal to the surface at the given test
                  points. May be NULL if needsTestNormal() returns true.
      \param[in]  trialNormals
                  Pointer to a Fortran-ordered 2D array of dimensions
                  (worldDimension(), pointCount) storing the components of
                  external unit vectors normal to the surface at the given test
                  points. May be NULL if needsTrialNormal() returns true.
      \param[out] result
                  Pointer to a preallocated Fortran-ordered 3D array of
                  dimensions (codomainDimension(), domainDimension(),
                  pointCount), in which will be stored the components of the
                  kernel function evaluated at the given point pairs. */
    void evaluate(int pointCount,
                  const CoordinateType* testPoints,
                  const CoordinateType* trialPoints,
                  const CoordinateType* testNormals,
                  const CoordinateType* trialNormals,
                  ValueType* result) const {
        asImp().evaluate(pointCount, testPoints, trialPoints,
                         testNormals, trialNormals, result);
    }

private:
    /**  Barton-Nackman trick */
    KernelImp& asImp() {
        return static_cast<KernelImp&>(*this);
    }

    /**  Barton-Nackman trick */
    const KernelImp& asImp() const {
        return static_cast<const KernelImp&>(*this);
    }
};

} // namespace Fiber

#endif
