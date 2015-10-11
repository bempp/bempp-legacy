Quadrature
==========

At its core the assembly of boundary integral operators in BEM++ is based on
an efficient evaluation of integrals of the form

.. math::

    I = \int_{T_1}\int_{T_2}g(x,y)\phi(y)\overline{\psi(x)}ds(y)ds(x),

where :math:`g(x,y)` is a weakly singular kernel and :math:`T_1` and :math:`T_2` are
triangles. :math:`\psi(x)` and :math:`\phi(y)` are basis functions on :math:`T_1` and 
:math:`T_2`. We need to distinguish two cases:

    * :math:`T_1` and :math:`T_2` are fully disjoint,
    * :math:`T_1` and :math:`T_2` share a vertex, edge or are identical.

In the first case a fully regular Gauss quadrature rule on triangles can be used, that is we 
approximate the integral as

.. math::

    I\approx \sum_{i=1}^{N_1}\sum_{j=1}^{N_2}g(x_i,y_j)\phi(y_j)\overline{\psi(x_i)}\omega^{(1)}_i\omega^{(2)}_j.

In the second case the singularity needs to be taken into account. This is done using Duffy type transformations.
Details of the singular integration can be found in the book `Boundary Element Methods <http://www.springer.com/us/book/9783540680925>`_
by Sauter and Schwab.

How many integration points are used?
-------------------------------------

The cost of the assembly of boundary element matrices is dominated by the number of kernel evaluations
for each pair of triangles :math:`(T_1,T_2)` in the approximation of the corresponding Galerkin integral.

For disjoint triangles a standard tensor Gauss rule based on symmetric Gauss points over triangles is used. The following table
shows the order, the number of points used on a triangle and the total number of kernel evaluations over the pair
:math:`(T_1, T_2)`.

====== ======================== ==================================
Order  Gauss points in triangle Total number of kernel evaluations 
====== ======================== ==================================
1      1                        1
2      3                        9
3      4                        16
4      6                        36
5      7                        49
6      12                       144
7      13                       169
8      16                       256
9      19                       361 
10     25                       625
11     27                       729
12     33                       1089
13     37                       1369
14     42                       1764
15     48                       2304
16     52                       2704
17     61                       3721
18     70                       4900
19     73                       5329
20     79                       6241
====== ======================== ==================================

If the triangles :math:`T_1` and :math:`T_2` are not disjoint the total number of kernel evaluations
depends on how the two triangles intersect each other. Let :math:`O` be the singular integration order.
Then the total number of kernel evaluations :math:`E_K` is given as follows.

    * The triangles are adjacent at a single vertex: :math:`E_K = 2 \left\lfloor (O+2)/2\right\rfloor^4`
    * The triangles are adjacent at an edge: :math:`E_K = 5 \left\lfloor (O+2)/2\right\rfloor^4`
    * The triangles coincide (:math:`T_1==T_2`) :math:`E_K = 6 \left\lfloor (O+2)/2\right\rfloor^4`

The default singular quadrature order in BEM++ is 6. Hence, the total number of kernel evaluations for an
element pair is 512 in the case of adjacent vertices, 1280 in the case of adjacent edges and 1536 in the
case of coincident triangles. This is significantly more than in the case of disjoint triangles. However,
the number of non-disjoint triangle pairs depends only linearly on the mesh size and typically only take
a small fraction of the overall computing time.

Controlling the quadrature order in BEM++
-----------------------------------------

BEM++ provides fine-grained control over the quadrature order. The simplest way to change the quadrature
order is to use the global parameter object ``bempp.api.global_parameters``. For example, to change the
singular quadrature order to :math:`5` use the command

::

    bempp.api.global_parameters.quadrature.double_singular = 5

The quadrature order for approximating the regular integral over pairs of disjoint triangles :math:`T_1` and
:math:`T_2` is dependent on the normalized distance between them. The normalized distance is defined as
the distance of the centers of the triangles divided by the maximum of the areas of :math:`T_1` and :math:`T_2`.

Three integration zones are defined, `near`, `medium` and `far`. The following parameters are set as default.

====== =================== ============ ============
Zone   Normalized distance Double Order Single Order
====== =================== ============ ============
near   2                   4            4
medium 4                   3            3 
far    infinity            2            2
====== =================== ============ ============

The given normalized distance is the maximum distance until which a quadrature rule is valid. The single order
is the quadrature order used for the evaluation of potentials, where the normalized distance is the relative
distance of an evaluation point to the corresponding surface element. For the discretisation of grid
functions the `far` single order is used. The values given here are relatively sane defaults. However, we
note that good values for the quadrature orders depend on the geometry and the accuracy requirements and
may be very different from the default values stated here.

To change for example the values of the medium distance quadrature orders the following commands can be used.

::

    bempp.api.global_parameters.quadrature.medium.max_rel_dist = 5
    bempp.api.global_parameters.quadrature.medium.single_order = 1
    bempp.api.global_parameters.quadrature.medium.double_order = 1


