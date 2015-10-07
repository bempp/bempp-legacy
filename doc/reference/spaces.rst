Function spaces
===============

Similarly to finite element codes function spaces form a central part of
BEM++. To initialize a function space all we need is a grid object.

.. code:: python

    import bempp.api
    grid = bempp.api.shapes.regular_sphere(3)

Let us now create a space of piecewise constant functions on the
elements. This is the standard low-order space to represent Neumann data
(that is normal derivatives) in boundary element methods.

.. code:: python

    space = bempp.api.function_space(grid, "DP", 0)

The first parameter of the ``function_space`` function is always the
grid. The second gives the type of space, in this case "DP" for
*D*\ iscontinuous-\ *P*\ olynomial, and the third parameter is the order
of the space (0 for piecewise constant).

We can now query the degrees of freedom of the space.

.. code:: python

    number_of_global_dofs = space.global_dof_count
    print("Number of global dofs: {0}".format(number_of_global_dofs))


.. parsed-literal::

    Number of global dofs: 512


For this space it is identical to the number of elements on the mesh.
This is not necessarily always the case.

Let us now create a space of continuous, piecewise linear functions.
This is the standard low-order space to represent Dirichlet data.

.. code:: python

    space = bempp.api.function_space(grid, "P", 1)
    number_of_global_dofs = space.global_dof_count
    print("Number of global dofs: {0}".format(number_of_global_dofs))


.. parsed-literal::

    Number of global dofs: 258


The number of global dofs is now identical to the number of vertices in
the grid.

Types of spaces
---------------

BEM++ supports the following types of spaces. The identifier for the
``function_space`` function is given in brackets.

-  Discontinuous polynomial spaces (DP). These are spaces of functions
   that are polynomial on each element but discontinuous across
   elements. The maximum order is 10.

-  Polynomial spaces (P). These are spaces of functions that are
   polynomial on each element and continuous across elements. The
   minimum order is zero. The maximum order is 10.

-  Polynomial spaces on barycentric grids (B-P). These are the same
   spaces as with the "P" identifier and the same number of degrees of
   freedom. But the underlying grid is a barycentric refinement of the
   original grid passed to the function (the barycentric refinement is
   created internally). This is needed in operator preconditioning.
   Currently, only ``order == 1`` is supported.
-  Discontinuous polynomial spaces on barycentric grids (B-DP). As
   above, but corresponding to discontinuous polynomial spaces on the
   original grid. Currently, only ``order == 1`` is supported.
-  Dual spaces of constant functions (DUAL). This is a space of constant
   functions defined on a dual grid (the dual grid is created internally
   from the grid object). These spaces form a stable dual pairing
   together with continuous, piecewise linear functions and are needed
   for certain opposite order preconditioners.
-  Raviart-Thomas Vector Space ("RT"). These are spaces of
   Raviart-Thomas basis functions. They are needed for integral
   operators in electromagnetic scattering. Currently, only low-order
   (``order == 0``) Raviart-Thomas spaces are supported.

For most scalar applications piecewise constant and continuous,
piecewise linear spaces are sufficient. For the electric field and
magnetic field operators spaces of Raviart-Thomas functions are
required. The barycentric and dual spaces are for the assembly of
certain types of preconditioned operators.

Local and global dofs
---------------------

An important concept for spaces are global and local degrees of freedom.
Global degrees of freedom are the actual dofs of the space associated
with the basis functions, while local dofs are the coefficients of the
restrictions of the basis functions to individual elements. Consider for
example a space of continuous, piecewise linear basis functions. Each
global dof is associated with a vertex. The corresponding basis function
lives on all elements that are adjacent to this vertex. Sometimes it is
useful to find out to what global dofs the basis functions on an element
contribute. This is shown in the following example.

.. code:: python

    space = bempp.api.function_space(grid, "P", 1)
    elements = list(grid.leaf_view.entity_iterator(0))
    element = elements[0]
    global_dofs, weights = space.get_global_dofs(element, dof_weights=True)
    print("Map from local to global dofs on element: {0}".format(global_dofs))
    print("Dof weights: {0} ".format(weights))


.. parsed-literal::

    Local2GlobalDofMap on element: [2, 66, 68]
    Dof weights: [1.0, 1.0, 1.0] 


The map from local to global dofs on the element denotes which local
basis function is assocated with which global basis functions. In this
example, the element has three local basis functions (associated with
the three vertices of the element), and the values of ``global_dofs``
are the indices of the global basis functions that they map to. The
array ``weights`` returns scalar multipliers that are the weights with
which the local basis function contributes to the global basis function.
Hence, to obtain the local coefficient of a local basis function the
corresponding global dof coefficient needs to be multiplied with the
``weight``. For most basis types the weight is always 1. But for
example, for Raviart-Thomas basis functions it can differ. Weights are
only returned if the parameter ``dof_weights=True`` is set.

Class descriptions
------------------

Function spaces are instantiated by the function :class:`bempp.api.function_space`.
They return :class:`bempp.api.space.Space` objects, which contain the necessary information that define
a space over a grid.

.. autofunction:: bempp.api.function_space

.. autoclass:: bempp.api.space.Space
    :members:

