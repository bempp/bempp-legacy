Grid Functions
==============

In BEM++ data on a given grid is represented as a ``GridFunction``
object. A ``GridFunction`` consists of a set of basis function
coefficients and a corresponding ``Space`` object. In the following we
will discuss the different ways of creating a grid function in BEM++.

Initializing a grid function from a Python callable
---------------------------------------------------

Often, in applications we are given an analytic expression for boundary
data. This is for example the case in many acoustic scattering problems,
where a typical scenario is that the incoming data is a plane wave or a
sound source. The following code defines a wave travelling with unit
wavenumber along the x-axis in the positive direction.

.. code:: python

    import numpy as np
    def fun(x, normal, domain_index, result):    
        result[0] = np.exp(1j * x[0])

A valid Python callable always has the same interface as shown above.
The first argument ``x`` is the coordinate of an evaluation point. The
second argument ``normal`` is the normal direction at the evaluation
point. The third one is the ``domain_index``. This corresponds to the
physical id in Gmsh and can be used to assign different boundary data to
different parts of the grid. The last argument ``result`` is the
variable that stores the value of the callable. It is a numpy array with
as many components as the basis functions of the underlying space have.
If a scalar space is given then ``result`` only has one component.

In order to discretise this callable we need to define a suitable space
object. Below we define a space of continuous, piecewise linear
functions on a spherical grid.

.. code:: python

    import bempp.api
    grid = bempp.api.shapes.regular_sphere(5)
    space = bempp.api.function_space(grid, "DP", 1)

The next command now discretizes the Python callable by projecting it
onto the space.

.. code:: python

    grid_fun = bempp.api.GridFunction(space, fun=fun)

There is quite a lot going on now as shown by the output logging
messages. Before we describe in detail what is happening we want to
visualize the ``grid_fun`` object. This can be done with the command

::

    grid_fun.plot()

This command opens ``Gmsh`` externally as a viewer to show the
``GridFunction`` object. The results should look like the following.

.. image:: figures/grid_fun.png

By default the real part of ``grid_fun`` is plotted. There are more
advanced functions to control this behavior.

Let us take a closer look at what happens in the initialization of this
GridFunction. Denote the global basis functions of the space by
:math:`\Psi_j`, :math:`j=1,\dots,N`. The computation of the grid
function consists of two steps.

1. Compute the projection coefficients
   :math:`p_j = \int_{\Gamma}\overline{\Psi_j}(\xi)f(\xi)d\xi`, where
   :math:`f` is the analytic function to be converted into a grid
   function and :math:`\Gamma` is the surface defined by the grid.
2. Compute the basis coefficients :math:`c_j` from :math:`Mc=p`, where
   :math:`M` is the mass matrix defined by
   :math:`M_{ij} = \int_{\Gamma}\overline{\Psi_i}(\xi)\Psi_j(\xi)d\xi`.

This is an orthogonal :math:`L^2(\Gamma)`-projection onto the basis
given by the :math:`\Psi_j`.

Initializing a grid function from coefficients or projections
-------------------------------------------------------------

Instead of an analytic expression we can initialize a ``GridFunction``
object also from a vector ``c`` of coefficients or a vector ``p`` of
projections. This can be done as follows.

::

    grid_fun = GridFunction(space, coefficients=c)
    grid_fun = GridFunction(space, projections=p, dual_space=dual)

The argument dual\_space gives the space with which the projection
coefficients were computed. The parameter is optional and if it is not
given then ``space == dual_space`` is assumed (i.e.
:math:`L^2(\Gamma)`-projection).

Function and class reference
----------------------------

.. autoclass:: bempp.api.GridFunction
    :members:

