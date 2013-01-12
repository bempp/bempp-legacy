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

/** \defgroup maxwell_3d Maxwell equations in 3D
    \ingroup operators

    This submodule contains classes relevant to the implementation of boundary
    operators and potential operators related to the time-harmonic Maxwell
    equations in 3D. Inside a region with a uniform scalar permittivity and
    permeability, these equations can be reduced to the form

    \f[
        \boldsymbol\nabla \times \boldsymbol\nabla \times
        \boldsymbol u - k^2 \boldsymbol u = 0,
    \f]

    where \f$\boldsymbol u\f$ stands either for the electric or the magnetic
    field and the number \f$k \neq 0\f$ is referred to as the <em>wave
    number</em>.

    The treatment of Maxwell equations in BEM++ closely follows that of A. Buffa
    and R. Hiptmair, "Galerkin Boundary Element Methods for Electromagnetic
    Scattering", in "Computational Methods in Wave Propagation", M. Ainsworth,
    ed., Springer, New York 2003, pp. 85-126. In the following we cite the
    principal definitions relevant to the implementation of BEM in the Maxwell
    case.

    The starting point is the Stratton-Chu theorem (Theorem 6 in Buffa and
    Hiptmair), which plays a part analogous to the Green's representation
    theorem for the Laplace equation. It states that any solution \f$u\f$ of the
    Maxwell equations in a bounded domain \f$\Omega\f$ with boundary
    \f$\Gamma\f$ has the representation

    \f[ \boldsymbol u =
    \boldsymbol \Psi_{\mathrm{SL}}(\gamma_{\mathrm{N,int}}\boldsymbol u) +
    \boldsymbol \Psi_{\mathrm{DL}}(\gamma_{\mathrm{D,int}}\boldsymbol u).
    \f]

    (all the symbols will be defined shortly). Furthermore, any solution \f$u\f$
    that satisfies the Maxwell equations in a domain \f$\Omega^{\mathrm c}\f$
    extending to infinity and meets the Silver-Mueller radiation conditions at
    infinity has the representation

    \f[ \boldsymbol u =
    -\boldsymbol \Psi_{\mathrm{SL}}(\gamma_{\mathrm{N,ext}}\boldsymbol u) -
    \boldsymbol \Psi_{\mathrm{DL}}(\gamma_{\mathrm{D,ext}}\boldsymbol u).
    \f]

    The action of the <em>single-layer potential operator</em> is defined by

    \f[
    (\boldsymbol \Psi_{\mathrm{SL}} \boldsymbol v)(\boldsymbol x) \equiv
    \mathrm i k \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    \boldsymbol v(\boldsymbol y) \mathrm\Gamma(\boldsymbol y) -
    \frac{1}{\mathrm i k}
    \boldsymbol\nabla_{\boldsymbol x}
    \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    (\boldsymbol \nabla_\Gamma \cdot \boldsymbol v)(\boldsymbol y)
    \mathrm\Gamma(\boldsymbol y),
    \f]

    where

    \f[
    G(\boldsymbol x, \boldsymbol y) \equiv
    \frac{\exp(\mathrm i k \lvert \boldsymbol x - \boldsymbol y \rvert)}
    {4\pi \lvert \boldsymbol x - \boldsymbol y \rvert}
    \f]

    is the Green's function of the Helmholtz equation. The action of the
    <em>double-layer potential operator</em>, in turn, is defined by
    \f[
    (\Psi_{\mathrm{DL}} \boldsymbol v)(\boldsymbol x) \equiv
    \boldsymbol\nabla_{\boldsymbol x} \times
    \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    \boldsymbol v(\boldsymbol y) \mathrm\Gamma(\boldsymbol y).
    \f]

    The <em>interior Dirichlet trace</em> \f$\gamma_{\mathrm{D,int}}\boldsymbol
    u\f$ at a point \f$\boldsymbol x \in \Gamma\f$ is defined as

    \f[(\gamma_{\mathrm{D,int}}\boldsymbol u)(\boldsymbol x) \equiv
    \boldsymbol u\rvert_{\Gamma,\mathrm{int}}(\boldsymbol x) \times
    \boldsymbol n(\boldsymbol x),\f]

    where \f$\boldsymbol n\f$
    is the outward unit vector normal to \f$\Gamma\f$ at \f$\boldsymbol x\f$ and
    \f$\boldsymbol u\rvert_{\Gamma,\mathrm{int}}(\boldsymbol x)\f$
    is the limit of \f$\boldsymbol u(\boldsymbol y)\f$ as
    \f$\boldsymbol y\f$ approaches \f$\boldsymbol x\f$
    from within \f$\Omega\f$.
    The <em>interior Neumann trace</em> \f$\gamma_{\mathrm{N,int}}\boldsymbol
    u\f$ at \f$\boldsymbol x \in \Gamma\f$ is defined as

    \f[(\gamma_{\mathrm{N,int}}\boldsymbol u)(\boldsymbol x) \equiv
    \frac{1}{\mathrm i k}
    (\boldsymbol \nabla \times \boldsymbol u)\rvert_{\Gamma,\mathrm{int}}
    (\boldsymbol x) \times
    \boldsymbol n(\boldsymbol x).\f]

    Note that compared to Buffa and Hiptmair, we include an additional
    \f$\mathrm i\f$ factor in the denominator of the Neumann trace, as this
    makes both the Neumann trace and the potential operators real-valued in the
    special case of \f$k\f$ purely imaginary (and \f$\boldsymbol u\f$ purely
    real).

    The <em>exterior</em> Dirichlet and Neumann traces are defined analogously,
    with the relevant quantities assumed to approach the point \f$\boldsymbol
    x\f$ from within the complement of \f$\Omega\f$.

    The <em>single-layer boundary operator</em>
    \f$\boldsymbol S\f$ is defined as the average of the interior and exterior
    Dirichlet traces of the single-layer potential operator:

    \f[
    \boldsymbol S =
    \frac12
    (\gamma_{\mathrm{D,int}} \boldsymbol \Psi_{\mathrm{SL}} +
    \gamma_{\mathrm{D,ext}} \boldsymbol \Psi_{\mathrm{SL}}).
    \f]

    The <em>double-layer boundary operator</em> \f$\boldsymbol C\f$ is defined,
    analogously, as the average of the interior and exterior Dirichlet traces of
    the double-layer potential operator.

    To exploit the symmetry between the electric and magnetic field in Maxwell
    equations, it is advantageous to use the following antisymmetric
    pseudo-inner product in the Galerkin discretisation of the above boundary
    operators:

    \f[
    \langle \boldsymbol u, \boldsymbol v \rangle_{\boldsymbol\tau,\Gamma} \equiv
    \int_\Gamma \boldsymbol u^* \cdot (\boldsymbol v \times \boldsymbol n)
    \mathrm d\Gamma,
    \f]

    where \f$*\f$ denotes complex conjugation. The weak forms of the single- and
    double-layer boundary operators with respect to this pseudo-inner product
    are given by

    \f[
    \langle \boldsymbol u,
    \boldsymbol S \boldsymbol v\rangle_{\boldsymbol\tau,\Gamma} =
    \int_\Gamma \int_\Gamma
    G(\boldsymbol x, \boldsymbol y)
    \biggl[-\mathrm i k
    \boldsymbol u^*(\boldsymbol x) \cdot \boldsymbol v(\boldsymbol y)
    -\frac{1}{\mathrm i k}
    (\boldsymbol \nabla_\Gamma \cdot \boldsymbol u^*)(\boldsymbol x)
    (\boldsymbol \nabla_\Gamma \cdot \boldsymbol v)(\boldsymbol y)
    \biggr]
    \mathrm{d}\Gamma(\boldsymbol x)\,\mathrm{d}\Gamma(\boldsymbol y)
    \f]

    and

    \f[
    \langle \boldsymbol u,
    \boldsymbol C \boldsymbol v\rangle_{\boldsymbol\tau,\Gamma} =
    \int_\Gamma \int_\Gamma
    \boldsymbol\nabla_{\boldsymbol x}G(\boldsymbol x, \boldsymbol y) \cdot
    [\boldsymbol u^*(\boldsymbol x) \times \boldsymbol v(\boldsymbol y)]
    \mathrm{d}\Gamma(\boldsymbol x)\,\mathrm{d}\Gamma(\boldsymbol y).
    \f]
 */

////////////////////////////////////////////////////////////////////////////////

/** \defgroup space Space
 *
 *  This module contains classes representing function spaces.
 */

#endif // bempp_doxygen_main_hpp
