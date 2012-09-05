The ``bempp.core`` module
=========================

Introduction
------------

This module contains SWIG-generated wrappers of C++ classes from BEM++.

Many of the classes occur in several variants corresponding to particular
instantiations of C++ class templates depending on parameters such as
``BasisFunctionType`` and ``ResultType`` (see the documentation of the C++
interface for more information on the meaning of these types). These variants
are distinguished by their names. For example, the group of classes denoted
symbolically as ``BoundaryOperator_BasisFunctionType_ResultType`` has six
concrete representatives:

- ``BoundaryOperator_float32_float32``,
- ``BoundaryOperator_float32_complex64``,
- ``BoundaryOperator_complex64_complex64``,
- ``BoundaryOperator_float64_float64``,
- ``BoundaryOperator_float64_complex128`` and
- ``BoundaryOperator_complex128_complex128``.

Here, ``float32`` and ``float64`` stand for single- and
double-precision floating-point numbers, whereas ``complex64`` and
``complex128`` denote single- and double-precision complex numbers.

It would be tedious and error-prone to have to specify these numeric types
explicitly whenever a new object is constructed. For this reason, the
``bempp.lib`` module provides ``create...`` helper functions that construct new
instances of classes from ``bempp.core`` with automatically determined values of
``BasisFunctionType`` and ``ResultType``. These "non-member constructors" should
be used in preference to the usual constructors.

The documentation of this module is in preliminary stage. The ``Grid``-related
classes have been almost fully documented; for the remaining classes only member
signatures are provided. Please refer to the documentation of the C++ interface
for more complete information.

Reference
---------

.. automodule:: bempp.core

