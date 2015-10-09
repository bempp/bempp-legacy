
Operator concepts
=================

Boundary operators
------------------

The principle operator concept in BEM++ is that of a boundary operator.
A boundary operator

.. math::


   A: \mathcal{D}\rightarrow \mathcal{R},

is a mapping from a *domain* space :math:`\mathcal{D}` into a *range*
space :math:`\mathcal{R}`, where both :math:`\mathcal{D}` and
:math:`\mathcal{R}` are defined on a given surface grid. BEM++ does not
directly work with the boundary operator :math:`A` itself but with a
weak form

.. math::


   a(u, v) := \int_{\Gamma} [Au](\mu)\overline{v(\mu)}d\mu,\quad u\in\mathcal{D},~v\in\mathcal{V}

where :math:`\mathcal{V}` is the dual space to the range space
:math:`\mathcal{R}` (in BEM++ we use the keyword *dual\_to\_range* for
the space :math:`\mathcal{V}`).

Boundary operators are defined in the subpackage
``bempp.api.operators.boundary``. Apart from Maxwell operators all
boundary operators take three space arguments: ``domain``, ``range``,
and ``dual_to_range``, which correspond to the spaces
:math:`\mathcal{D}`, :math:`\mathcal{R}` and :math:`\mathcal{V}`. The
following code snippet defines the Laplace single-layer boundary
operator on a space of piecewise constant functions. For simplicity, we
choose all three space arguments to be identical.

.. code:: python

    import bempp.api
    grid = bempp.api.shapes.regular_sphere(3)
    space = bempp.api.function_space(grid, "DP", 0)
    slp = bempp.api.operators.boundary.laplace.single_layer(space, space, space)


.. parsed-literal::

    INFO:BEMPP:Found Dolfin. FEM/BEM coupling with FEniCS enabled.


It is important to note that the above code only sets up data
structures. No discretisation of an operator is happening.

A complete algebra for operators is implemented. We can add operators,
multiply them with scalars, and also multiply operators. Hence, the
following operations are all valid.

.. code:: python

    scaled_operator = 1.5 * slp
    sum_operator = slp + slp
    squared_operator = slp * slp

Particularly, interesting is the last step. Assume that the matrix
:math:`S` is the Galerkin discretisation of the ``slp`` operator with
the given space of piecewise constant functions. Then the discretisation
of ``squared_operator`` is computed as :math:`SM^{-1}S`, where :math:`M`
is the mass matrix of inner products between functions in the
``dual_to_range`` space and the ``range`` space. This is done
automatically in BEM++ so that the user does not have to deal with the
correct mass matrix operations manually. It allows a complete product
algebra based on Galerkin discretisations. Note that if the mass matrix
:math:`M` is not square, BEM++ internally solves a least-squares problem
using a normal equation approach.

Operators can also be multiplied with grid functions as shown in the
following.

.. code:: python

    # Create a grid function with unit coefficients.
    import numpy as np
    number_of_global_dofs = space.global_dof_count
    coeffs = np.ones(number_of_global_dofs)
    grid_fun = bempp.api.GridFunction(space, coefficients=coeffs)
    
    # Now apply the operator to the grid function
    
    result_fun = slp * grid_fun

There is quite a lot going on internally for this operation as shown by
the logging messages. First of all, in order to apply the operator
``slp`` to the grid function ``grid_fun`` BEM++ needs to assemble the
operator. This is the first step, where an actual discretisation is
computed. However, assembling the operator is not enough. To compute the
coefficients :math:`c_{new}` of the grid funtion ``result_fun`` a mass
matrix needs to be assembled since the coefficients of the result are
computed as

.. math::


   c_{new} = M^{-1}Sc,

where :math:`c` is the vector of coefficients of ``grid_fun``. Again,
all this is handled automatically by BEM++.

Discrete boundary operators
---------------------------

Quite often it is necessary to have more direct access to
discretisations of boundary operators. For this a boundary operator
provides two methods ``weak_form`` and ``strong_form``. Given the
variational form :math:`a(u,v)` as defined above the discretisation of
this variational form is the matrix :math:`A_h` defined by

.. math::


   [A_h]_{i,j} = a(\psi_i, \phi_j),

where the :math:`\psi_i` are the basis functions of the test space
:math:`\mathcal{V}` and the :math:`\phi_j` are the basis functions of
the domain space :math:`\mathcal{D}`.

Discrete boundary operators give access to the discretized matrix by
providing routines to perform matrix-vector products and query the
underlying matrix. The following gives an example.

.. code:: python

    slp_discrete = slp.weak_form()
    print("Shape of the matrix: {0}".format(slp_discrete.shape))
    print("Type of the matrix: {0}".format(slp_discrete.dtype))
    
    x = np.random.rand(slp_discrete.shape[1])
    y = slp_discrete * x


.. parsed-literal::

    Shape of the matrix: (512, 512)
    Type of the matrix: float64


Discrete boundary operators implement the ``LinearOperator`` protocol
provided by recent ``SciPy`` versions. This means that the standard
operations such as multiplications with scalars, addition, and
multiplication are available.

It is possible to convert a discrete boundary operator into a standard
``NumPy`` array. However, this is not recommended for larger problems.
Depending on the discretisation mode a discrete operator may store a
matrix only implicitly and not as a dense array. This then needs to be
converted to a dense array by multiplication with an identity matrix.
The following command turns a discrete boundary operator into a
``NumPy`` array.

.. code:: python

    slp_mat = bempp.api.as_matrix(slp_discrete)
    print(slp_mat)


.. parsed-literal::

    [[  6.01598524e-04   1.59955449e-04   1.60210104e-04 ...,   3.15819132e-05
        3.35583598e-05   3.43634489e-05]
     [  1.59811402e-04   6.20059370e-04   1.21665497e-04 ...,   3.14103888e-05
        3.27991479e-05   3.39030174e-05]
     [  1.60210104e-04   1.21313696e-04   6.20059370e-04 ...,   3.33189379e-05
        3.55234393e-05   3.65203446e-05]
     ..., 
     [  3.15574627e-05   3.14105754e-05   3.33145031e-05 ...,   1.68183155e-03
        4.05004969e-04   7.46497936e-04]
     [  3.35100884e-05   3.27972380e-05   3.54902004e-05 ...,   4.05005071e-04
        1.68183155e-03   7.46490617e-04]
     [  3.43687317e-05   3.39000037e-05   3.65254704e-05 ...,   7.46498073e-04
        7.46490752e-04   1.79722275e-03]]


Above we have used the method ``weak_form`` to obtain the discretisation
of a boundary operator. There is another method called ``strong_form``.
While the method ``weak_form()`` returns a discretisation of the
variational form :math:`a(u,v)` the method ``strong_form`` returns a
discretisation of the action of the original operator :math:`A` into the
range space. Hence, the operators returned by the two methods are as
follows:

-  weak\_form: :math:`[A_h]_{i,j} = a(\psi_i, \phi_j)`
-  strong\_form: :math:`M^{-1}A_h`

Here, :math:`M` is the mass matrix between the ``range`` space and the
``dual_to_range`` space. We note that :math:`M^{-1}` is never formed
explicitly but the corresponding system is solved internally by sparse
LU decomposition, where the factorisation is only done once and then
stored.

Potential operators
-------------------

A potential operator, as opposed to a boundary operator in BEM++, maps
from a given space over the boundary grid :math:`\Gamma` into a set of
external evaluation points :math:`x_j\not\in\Gamma`. Let us demonstrate
at a simple example how to evaluate the Laplace single-layer potential
operator at certain points away from the boundary.

.. code:: python

    # Define two evaluation points
    evaluation_points = np.array([[2, 3],
                                  [1, 0],
                                  [4, 5]])
    
    # Create the Laplace single-layer potential operator
    slp_pot = bempp.api.operators.potential.laplace.single_layer(space, evaluation_points)
    potential_values = slp_pot * grid_fun
    print(potential_values)


.. parsed-literal::

    [[ 0.21547098  0.16933975]]


As opposed to boundary operators potentials are assembled immediately
when they are instantiated. Potential operators implement a simple
algebra. Hence, they allow multiplication with scalars and addition with
other potentials. To apply a potential to a given surface density it can
be multiplied with a grid function as shown above. The result is an
array of potential values, in which each column consist of the
components of the potential at a given evaluation point. In this case
the potential is scalar. Hence, each column has only one entry.

Function and class reference
----------------------------

.. autoclass:: bempp.api.assembly.BoundaryOperator
    :members:

.. autoclass:: bempp.api.assembly.ZeroBoundaryOperator

.. autoclass:: bempp.api.assembly.ElementaryBoundaryOperator

.. autoclass:: bempp.api.assembly.LocalBoundaryOperator

.. autoclass:: bempp.api.assembly.GeneralNonlocalDiscreteBoundaryOperator
    :members:

.. autoclass:: bempp.api.assembly.DenseDiscreteBoundaryOperator
    :members:

.. autoclass:: bempp.api.assembly.SparseDiscreteBoundaryOperator
    :members:

.. autoclass:: bempp.api.assembly.InverseSparseDiscreteBoundaryOperator

.. autoclass:: bempp.api.assembly.ZeroDiscreteBoundaryOperator

.. autofunction:: bempp.api.assembly.as_matrix

.. autoclass:: bempp.api.assembly.PotentialOperator
    :members:





