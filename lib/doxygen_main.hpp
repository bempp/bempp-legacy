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

////////////////////////////////////////////////////////////////////////////////

/** \defgroup assembly Assembly
 *
 *  This is the most important module of BEM++. It contains classes
 *  representing boundary and potential operators, classes representing
 *  functions defined on surfaces and in volumes, and classes related to the
 *  assembly of weak forms of boundary operators.
 */

/** \defgroup abstract_boundary_operators Abstract boundary operators
 *  \ingroup assembly
 *
 *  This submodule contains the AbstractBoundaryOperator class and its
 *  descendants. Note: subclasses of AbstractBoundaryOperator related to
 *  particular differential equations are included in the \ref operators module.
 */

/** \defgroup boundary_operators Boundary operators
 *  \ingroup assembly
 *
 *  This submodule contains classes that encapsulate one or more abstract
 *  boundary operators together with their weak forms.
 */

/** \defgroup discrete_boundary_operators Discrete boundary operators
 *  \ingroup assembly
 *
 *  This submodule contains the DiscreteBoundaryOperator class and its
 *  descendants.
 */

/** \defgroup potential_operators Potential operators
 *  \ingroup assembly
 *
 *  This submodule contains the PotentialOperator class and its descendants.
 *  Note: subclasses of PotentialOperator related to particular differential
 *  equations are included in the \ref operators module.
 */

/** \defgroup assembly_functions Functions
 *  \ingroup assembly
 *
 *  This submodule contains classes representing functions defined on surfaces
 *  and in volumes.
 */

/** \defgroup weak_form_assembly Weak-form assembly
 *  \ingroup assembly
 *
 *  This submodule contains classes controlling weak-form assembly.
 */

/** \defgroup weak_form_assembly_internal Weak-form assembly (internal)
 *  \ingroup assembly
 *
 *  This submodule contains classes used internally by BEM++ during weak-form
 *  assembly.
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup common Common utilities
 *
 *  This module containts utility classes and functions used by other modules of
 *  BEM++.
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup fiber Fiber
 *
 *  This module contains classes implementing low-level quadrature and
 *  local (element-by-element) assembly.
 */

/** \defgroup quadrature Quadrature
 *  \ingroup fiber
 *
 *  This submodule contains classes responsible for evaluation of integrals.
 */

/** \defgroup weak_form_elements Weak-form elements
 *  \ingroup fiber
 *
 *  This submodule defines interfaces of constituent elements of weak forms of
 *  integral operators. It also contains default implementations of these
 *  interfaces.
 */

/** \defgroup functors
 *  \ingroup fiber
 *
 *  This submodule contains functor classes evaluating elements of the weak
 *  forms of specific operators at single points or pairs of points.
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup grid Grid
 *
 *  This module contains classes related to the management of boundary-element
 *  grids.
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup linalg Linear algebra
 *
 *  This module contains classes related to the solution of linear algebraic
 *  equations arising from discretizations of boundary integral equations.
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup operators Operators
 *
 *  This module contains classes implementing particular boundary operators and
 *  potential operators.
 */

/** \defgroup identity Identity operator
 *  \ingroup operators
 *
 *  This submodule contains the class representing the identity operator.
 */

/** \defgroup composite_boundary_operators Composite boundary operators
 *  \ingroup operators
 *
 *  This submodule contains classes representing composite abstract boundary
 *  operators, such as their sum.
 */

/** \defgroup laplace_3d Laplace equation in 3D
 *  \ingroup operators
 *
 *  This submodule contains classes implementing kernels, boundary operators and
 *  potential operators related to the Laplace equation in 3D,
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
 *  This submodule contains classes implementing kernels, boundary operators and
 *  potential operators related to the Helmholtz equation in 3D,
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
 *  This submodule contains classes implementing kernels, boundary operators and
 *  potential operators related to the modified Helmholtz equation in 3D,
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

////////////////////////////////////////////////////////////////////////////////

/** \defgroup space Space
 *
 *  This module contains classes representing function spaces.
 */

#endif // bempp_doxygen_main_hpp
