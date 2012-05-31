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

#ifndef bempp_doxygen_main_hpp
#define bempp_doxygen_main_hpp

// This file contains no code and therefore does not need to be included
// anywhere; its only purpose is to hold documentation that does not belong
// anywhere else.

/** \defgroup assembly Assembly
 *
 *  This is the most important module of BEM++. It contains classes related to
 *  the assembly of weak forms of linear operators.
 */

/** \defgroup common Common utilities
 *
 *  This module containts utility classes and functions used by other modules of
 *  BEM++.
 */

/** \defgroup fiber Fiber
 *
 *  This module contains classes implementing low-level quadrature and
 *  local (element-by-element) assembly.
 */

/** \defgroup linalg Linear algebra
 *
 *  This module contains classes related to the solution of linear algebraic
 *  equations arising from discretisations of boundary integral equations.
 */

/** \defgroup operators Operators
 *
 *  This module contains classes implementing particular boundary integral
 *  operators.
 */

/** \defgroup space Space
 *
 *  This module contains classes representing function spaces.
 */

/** \defgroup laplace_3d Laplace equation in 3D
 *  \ingroup operators
 *
 *  This module contains classes implementing kernels and operators related to the
 *  Laplace equation in 3D,
 *  \f[
 *      \biggl(\frac{\partial^2}{\partial x^2} +
 *      \frac{\partial^2}{\partial y^2} +
 *      \frac{\partial^2}{\partial z^2}\biggr)
 *      u(x, y, z) = 0.
 *  \f]
 */

/** \defgroup helmholtz_3d Helmholtz equation in 3D
 *  \ingroup operators
 *
 *  This module contains classes implementing operators related to the
 *  Helmholtz equation in 3D,
 *  \f[
 *      \biggl(\frac{\partial^2}{\partial x^2} +
 *      \frac{\partial^2}{\partial y^2} +
 *      \frac{\partial^2}{\partial z^2} +
 *      k^2\biggr)
 *      u(x, y, z) = 0.
 *  \f]
 *  The number \f$k\f$ is referred to as the <em>wave number</em>.
 *
 *  \note The term <em>wave number</em> refers to different physical quantities
 *  in the Helmholtz equation and in the
 *  \ref modified_helmholtz_3d "modified Helmholtz equation".
 */

/** \defgroup modified_helmholtz_3d Modified Helmholtz equation in 3D
 *  \ingroup operators
 *
 *  This module contains classes implementing kernels and operators related to the
 *  modified Helmholtz equation in 3D,
 *  \f[
 *      \biggl(\frac{\partial^2}{\partial x^2} +
 *      \frac{\partial^2}{\partial y^2} +
 *      \frac{\partial^2}{\partial z^2} -
 *      k^2\biggr)
 *      u(x, y, z) = 0.
 *  \f]
 *  The number \f$k\f$ is referred to as the <em>wave number</em>.
 *
 *  \note The term <em>wave number</em> refers to different physical quantities
 *  in the modified Helmholtz equation and in the
 *  \ref helmholtz_3d "standard Helmholtz equation".
 */

#endif // bempp_doxygen_main_hpp
