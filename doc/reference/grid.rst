Grids 
=====


Grids in BEM++ are implementeded using the `Dune-Grid library <http://www.dune-project.org>`_
and the Python classes model parts of the underlying Dune data structures.

The provides :class:`Grid` class defines the interfaces to iterate through grids and
query geometric information. It cannot be directly instantiated. BEM++ rather provides
factory functions and a :class:`GridFactory` object to construct grids.

.. autoclass:: bempp.api.grid.Grid
    :members:


