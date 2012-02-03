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

#ifndef fiber_geometry_hpp
#define fiber_geometry_hpp

namespace Fiber
{

/** \brief Geometry interface.

  This class template is used as a base class for all kernel implementations.
  It uses the Barton-Nackman trick to ensure conformity to the interface.

  A "geometry", as meant here, is an encapsulation of the mapping from
  reference (local) to physical (global) coordinates on an element.

  \tparam CoordinateType Type used to represent components of local and global
                         coordinates (e.g. float or double).
  \tparam GeometryImp    Type of an implementation of the geometry interface. */
template <typename CoordinateType, typename GeometryImp>
class Geometry
{
public:
    /** \brief Element dimension */
    int dimension() const {
        return asImp().dimension();
    }

    /** \brief Dimension of the space containing the element */
    int worldDimension() const {
        return dimension() + 1;
    }

    /** \brief Evaluate requested geometrical quantities at prescribed points.

      \param[in]  pointCount
                  Number of points.
      \param[in]  localCoords
                  Pointer to a Fortran-ordered 2D array of dimensions
                  (dimension(), pointCount) storing local coordinates of points
                  at which the requested geometrical quantities should be
                  evaluated.
      \param[out] globalCoords
                  Pointer to a preallocated Fortran-ordered 2D array of
                  dimensions (worldDimension(), pointCount), in which will be
                  stored the global coordinates of the given points.
      \param[out] normals
                  Pointer to a preallocated Fortran-ordered 2D array of
                  dimensions (worldDimension(), pointCount), in which will be
                  stored the components of unit vectors normal to the element at
                  the given points.
      \param[out] integrationElements
                  Pointer to a preallocated array of length pointCount, in
                  which will be stored the values of \$f\sqrt{|\det(J_g^T
                  J_g)|}\f$ at the given points, with \f$ J_g\f$ denoting the
                  Jacobian matrix of the local->global coordinate
                  transformation described by this geometry.
      \param[out] jacobiansTransposed
                  Pointer to a preallocated Fortran-ordered 3D array of
                  dimensions (dimension(), worldDimension(), pointCount), in
                  which will be stored the transposed Jacobian matrices of the
                  local->global coordinate transformation described by this
                  geometry at the given points.
      \param[out] jacobianInversesTransposed
                  Pointer to a preallocated Fortran-ordered 3D array of
                  dimensions (worldDimension(), dimension(), pointCount), in
                  which will be stored the pseudoinverses of the transposed
                  Jacobian matrices of the local->global coordinate
                  transformation described by this geometry at the given points.

      Each of the output arguments can be set to NULL, in which case the
      corresponding geometric quantity is not evaluated.
      */
    void getInformation(
            int pointCount,
            const CoordinateType* localCoords,
            CoordinateType* globalCoords,
            CoordinateType* normals,
            CoordinateType* integrationElements,
            CoordinateType* jacobiansTransposed,
            CoordinateType* jacobianInversesTransposed) {
        asImp().getInformation(
                    pointCount, localCoords,
                    globalCoords, normals, integrationElements,
                    jacobiansTransposed, jacobianInversesTransposed);
    }

private:
    /**  Barton-Nackman trick */
    GeometryImp& asImp() {
        return static_cast<GeometryImp&>(*this);
    }

    /**  Barton-Nackman trick */
    const GeometryImp& asImp() const {
        return static_cast<const GeometryImp&>(*this);
    }
};

}

#endif
