
Available operators
===================

In the following we describe all operators available in BEM++.

Scalar non-local operators (Laplace, modified Helmholtz and Helmholtz)
----------------------------------------------------------------------

Let :math:`g(x,y)` be a given Green's function. We define the
single-layer potential operator and double layer-potential operator as
follows:

.. math::


   [\mathcal{S}\phi](x) &= \int_{\Gamma}g(x,y)\phi(y)ds(y),~x\in\mathbb{R}^3\backslash\Gamma\nonumber\\
   [\mathcal{K}\phi](x) &= \int_{\Gamma}\frac{\partial g(x,y)}{\partial\nu(y)}\phi(y)ds(y),~x\in\mathbb{R}^3\backslash\Gamma

From this we derive the following boundary operators:

-  Single-Layer boundary operator:
   :math:`[S\phi](x) = \int_{\Gamma}g(x,y)\phi(y)ds(y),~x\in\Gamma`
-  Double-Layer boundary operator:
   :math:`[K\phi](x) = \int_{\Gamma}\frac{\partial g(x,y)}{\partial\nu(y)}\phi(y)ds(y),~x\in\Gamma`
-  Adjoint double layer boundary operator:
   :math:`[K'\phi](x) = \int_{\Gamma}\frac{\partial g(x,y)}{\partial\nu(x)}\phi(y)ds(y),~x\in\Gamma`
-  Hypersingular boundary operator:
   :math:`[D\phi](x) = -\frac{\partial}{\partial \nu(x)}\int_{\Gamma}\frac{\partial g(x,y)}{\partial\nu(y)}\phi(y)ds(y),~x\in\Gamma`


The actual implementation of boundary operators in BEM++ are based on
variational formulations. Hence, the implementation for example of the
single-layer boundary operator is given by

.. math::


   s(\phi, \psi) = \int_{\Gamma}\overline{\psi(x)}\int_{\Gamma}g(x,y)\phi(y)ds(y)ds(x).

For the hypersingular operator integration a slightly different
formulation based on integration by parts of this variational form is
used. Details can be found for example in the book `Numerical
Approximation Methods for Elliptic Boundary Value
Problems <http://www.numerik.math.tu-graz.ac.at/~steinbach/FEMBEMEng.htm>`_
by O. Steinbach.

The different types of boundary operators are dependent on the choice of
the Green's function :math:`g`. In BEM++ the definitions are as follows.

-  Laplace (:math:`-\Delta u = 0`):

   .. math:: g(x,y) = \frac{1}{4\pi |x-y|}

-  Modified Helmholtz (:math:`-\Delta u+\omega^2 u = 0`):

   .. math:: g(x,y) = \frac{e^{-\omega |x-y|}}{4\pi|x-y|}

-  Helmholtz (:math:`-\Delta u - k^2 u = 0`):

   .. math:: g(x,y) = \frac{e^{ik |x-y|}}{4\pi|x-y|}

The boundary operators are defined in the following modules:

-  Laplace: ``bempp.api.operators.boundary.laplace``
-  Helmholtz: ``bempp.api.operators.boundary.helmholtz``
-  Modified Helmholtz: ``bempp.api.operators.boundary.modified_helmholtz``

In each of these modules the names of the corrresponding functions are
``single_layer``, ``double_layer``, ``adjoint_double_layer`` and
``hypersingular``.

The packages for the corresponding potential operators are:

-  Laplace: ``bempp.api.operators.potential.laplace``
-  Helmholtz: ``bempp.api.operators.potential.helmholtz``
-  Modified Helmholtz:
   ``bempp.api.operators.potential.modified_helmholtz``

The names of the corresponding potentials in each module are
``single_layer`` and ``double_layer``.

For example, a Helmholtz hypersingular boundary operator is instantiated
as

::

    hyp = bempp.api.operators.boundary.helmholtz.hypersingular(domain, range, dual_to_range, k)

Here, :math:`k` is the wavenumber as defined above.

In several applications, in particular transmission problems it is useful to define the multitrace operator

.. math::

    A = \begin{bmatrix} -K & S\\ D & K'\end{bmatrix}.

This operator has the property that 

.. math::

    \left[\frac{1}{2}I + A\right]\begin{bmatrix}u \\ u_n \end{bmatrix} = \begin{bmatrix}u \\ u_n\end{bmatrix}

if :math:`u` and :math:`u_n` are the Dirichlet, respectively Neumann data, of an interior solution to the corresponding
PDE. The operator :math:`\frac{1}{2}I + A` is also called a Calderon projector.
From this property it follows directly that :math:`(2A)^2 = I` and therefore :math:`A` is self-regularizing.

Currently, the operator :math:`A` is implemented using linear, piecewise continuous basis functions for the Dirichlet
data and piecewise constant basis functions on a dual grid for the Neumann data. This choice of basis functions allows
to perform the product :math:`A * A` in a stable way.

The following defines a multitrace operator for a Helmholtz problem with wavenumber :math:`k=1` and forms the interior
Calderon projector :math:`\frac{1}{2}I+A`.

::

    grid = bempp.api.shapes.regular_sphere(3)
    A = bempp.api.operators.helmholtz.multitrace_operator(grid, 1)
    ident = bempp.api.operators.boundary.sparse.multitrace_identity(grid)
    calderon = .5 * ident + A

The assembly of the multitrace operator :math:`A` requires the assembly of the single-layer and double-layer boundary
operators using discontinuous piecewise linear basis functions on a barycentric refinement of the grid. This
barycentric refinement has six times the number of elements as the original grid. This is necessary to assemble
a basis for piecewise constant basis functions on the dual grid.



Scalar, sparse operators
------------------------

BEM++ implements two sparse operators acting on scalar spaces, the
identity operator :math:`I` and the Laplace-Beltrami operator
:math:`-\Delta_{\Gamma}`. They are given in their variational form as
follows.

-  Identity Operator:

   .. math::


      m(\phi, \psi) = \int_{\Gamma} \overline{\Psi(y)} \phi(y)dy

-  Laplace-Beltrami Operator:

   .. math::


      k(\phi, \psi) = \int_{\Gamma} \overline{\nabla_{\Gamma}\Psi(y)}\cdot \nabla_{\Gamma}\phi(y)dy

The corresponding functions are ``bempp.api.boundary.sparse.identity``
and ``bempp.api.boundary.sparse.laplace_beltrami``.

Maxwell operators
-----------------

BEM++ supports the solution of the time-harmonic Maxwell equation of the
form

.. math::


   \nabla\times\nabla\times u - k^2 u = 0.

We define the following two potentials:

-  Maxwell electric field potential operator:

   .. math::


      [\mathcal{E}\phi](x) = ik\int_{\Gamma}g(x,y)\phi(y)ds(y)-\frac{1}{ik}\nabla_x\int_{\Gamma}g(x,y)(\nabla_{\Gamma}\cdot\phi)(y)ds(y)

-  Maxwell magnetic field potential operator:

   .. math::


      [\mathcal{M}\phi](x) = \nabla_x\times\int_{\Gamma}g(x,y)\phi(y)ds(y)

   The corresponding functions are
   ``bempp.api.operators.potential.maxwell.electric_field`` and
   ``bempp.api.operators.potential.maxwell.magnetic_field``.

The definition of the electric field operator given above includes a
factor :math:`i` in the nominator. This is less common, but has
implementational advantages in BEM++.

Based on the electric and magnetic field potential operators we can
derive the corresponding boundary operators. We will not give details
but just state the variational formulations of the operators.

-  Maxwell electric field boundary operator:

   .. math::


      s(\phi, \psi) = \int_{\Gamma}\int_{\Gamma}g(x,y)\left[-ik\overline{\psi(x)}\cdot\phi(y)-\frac{1}{ik}\left(\nabla_{\Gamma}\cdot\overline{\psi}\right)(x)\left(\nabla_{\Gamma}\cdot\phi\right)(y)\right]ds(x)ds(y)

-  Maxwell magnetic field boundary operator:

   .. math::


      c(\phi, \psi) = \int_{\Gamma}\int_{\Gamma}\nabla_xg(x,y)\left[\overline{\psi(x)}\times \phi(y)\right]ds(x)ds(y)

   The corresponding BEM++ functions are
   ``bempp.operators.boundary.maxwell.electric_field`` and
   ``bempp.operators.boundary.maxwell.magnetic_field``.

The discretisation of Maxwell operators uses an antisymmetric
pseudo-inner product defined as

.. math::


   i(\phi, \psi) = \int_{\Gamma}\overline{\psi(y)}\left(\phi(y)\times \nu(y)\right)ds(y).

For Maxwell problems this should be used instead of the standard
identity operator. The corresponding function in BEM++ is
``bempp.operators.boundary.sparse.maxwell_identity``.

Far Field Operators
-------------------

BEM++ implements far field operators for Maxwell and Helmholtz problems.

Approximate DtN and NtD operators
---------------------------------

BEM++ implements the OSRC approximations to the Dirichlet-to-Neumann (DtN) and Neumann-to-Dirichlet (NtD) maps
(see e.g. X. Antoine, M. Darbas,  Generalized Combined Field Integral Equations for the Iterative Solution of the Three-Dimensional Helmholtz Equation, Mathematical Modelling and Numerical Analysis 41 (1) (2007), pp. 147-167). The corresponding functions are :class:`bempp.api.boundary.helmholtz.osrc_ntd` and :class:`bempp.api.boundary.helmholtz.osrc_dtn`.

Function and class reference (Laplace)
--------------------------------------
.. autofunction:: bempp.api.operators.boundary.laplace.single_layer
.. autofunction:: bempp.api.operators.boundary.laplace.double_layer
.. autofunction:: bempp.api.operators.boundary.laplace.adjoint_double_layer
.. autofunction:: bempp.api.operators.boundary.laplace.hypersingular
.. autofunction:: bempp.api.operators.boundary.laplace.multitrace_operator
.. autofunction:: bempp.api.operators.boundary.laplace.interior_calderon_projector
.. autofunction:: bempp.api.operators.boundary.laplace.exterior_calderon_projector
.. autofunction:: bempp.api.operators.potential.laplace.single_layer
.. autofunction:: bempp.api.operators.potential.laplace.double_layer

Function and class reference (Modified Helmholtz)
-------------------------------------------------
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.single_layer
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.double_layer
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.adjoint_double_layer
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.hypersingular
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.multitrace_operator
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.interior_calderon_projector
.. autofunction:: bempp.api.operators.boundary.modified_helmholtz.exterior_calderon_projector
.. autofunction:: bempp.api.operators.potential.modified_helmholtz.single_layer
.. autofunction:: bempp.api.operators.potential.modified_helmholtz.double_layer


Function and class reference (Helmholtz)
----------------------------------------

.. autofunction:: bempp.api.operators.boundary.helmholtz.single_layer
.. autofunction:: bempp.api.operators.boundary.helmholtz.double_layer
.. autofunction:: bempp.api.operators.boundary.helmholtz.adjoint_double_layer
.. autofunction:: bempp.api.operators.boundary.helmholtz.hypersingular

.. autofunction:: bempp.api.operators.boundary.helmholtz.multitrace_operator
.. autofunction:: bempp.api.operators.boundary.helmholtz.interior_calderon_projector
.. autofunction:: bempp.api.operators.boundary.helmholtz.exterior_calderon_projector
.. autofunction:: bempp.api.operators.boundary.helmholtz.osrc_dtn
.. autofunction:: bempp.api.operators.boundary.helmholtz.osrc_ntd

.. autofunction:: bempp.api.operators.potential.helmholtz.single_layer
.. autofunction:: bempp.api.operators.potential.helmholtz.double_layer

.. autofunction:: bempp.api.operators.far_field.helmholtz.single_layer
.. autofunction:: bempp.api.operators.far_field.helmholtz.double_layer

Function and class reference (Maxwell)
--------------------------------------

.. autofunction:: bempp.api.operators.boundary.maxwell.electric_field
.. autofunction:: bempp.api.operators.boundary.maxwell.magnetic_field
.. autofunction:: bempp.api.operators.potential.maxwell.electric_field
.. autofunction:: bempp.api.operators.potential.maxwell.magnetic_field
.. autofunction:: bempp.api.operators.far_field.maxwell.electric_field
.. autofunction:: bempp.api.operators.far_field.maxwell.magnetic_field

Function and class reference (sparse operators)
-----------------------------------------------
.. autofunction:: bempp.api.operators.boundary.sparse.identity
.. autofunction:: bempp.api.operators.boundary.sparse.maxwell_identity
.. autofunction:: bempp.api.operators.boundary.sparse.laplace_beltrami
.. autofunction:: bempp.api.operators.boundary.sparse.multitrace_identity


