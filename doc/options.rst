Changing the BEM++ options
==========================

BEM++ uses a global options structure in ``bempp.api.global_parameters`` to control
various aspects of the library. In addition several functions and classes have a
``parameters`` attribute that can take a custom parameters object. This is useful
if for a certain function one wants to override the global parameters. To create
a custom parameters object juse use the command

::

    parameters = bempp.api.common.global_parameters()

In the following we give a description of all global parameter options.

Description of all global parameters
------------------------------------

* ``bempp.api.global_parameters.assembly.boundary_operator_assembly_type``:
  Controls wheter boundary operators are assembled in `dense` mode or in `hmat` mode.
  The default is `hmat` to use H-Matrix assembly for boundary operators.
* ``bempp.api.global_parameters.assembly.potential_operator_assembly_type``:
  Controls wheter potential operators are assembled in `dense` mode or in `hmat` mode.
  The default is `hmat` to use H-Matrix assembly for potential operators. H-Matrix
  assembly for exterior potentials is almost always preferable.
* ``bempp.api.global_parameters.assembly.enable_singular_integral_caching``:
  If set to True (default) singular integrals are pre-calculated and cached
  before the regular integrals are calculated. This usually gives a small
  speed advantage and should not need to be modified.
* ``bempp.api.global_parameters.assembly.enable_interpolation_for_oscillatory_kernels``:
  If set to True (default) Helmholtz type kernels (including Maxwell) are evaluated
  using a piecewise Hermite interpolation. This is significantly faster than evaluating the exponentials
  directly. 
* ``bempp.api.global_parameters.assembly.interpolation_points_per_wavelength``:
  The number of interpolation points per wavelength to be used (default 5000).
* ``bempp.api.global_parameters.quadrature``: Modify the quadrature options
  (see also :doc:`quadrature`).
* ``bempp.api.global_parameters.hmat.admissibility``:
  The type of admissibility condition. The default is `weak`, which declares
  a block cluster admissible if the column and range cluster bounding boxes do not
  intersect. The classical H-Matrix admissibility condition is obtained by
  setting this parameter to `strong`.
* ``bempp.api.global_parameters.hmat.coarsening``:: If True (default) enable
  automatic coarsening of H-Matrices.
* ``bempp.api.global_parameters.hmat.coarsening_accuracy``: The accuracy
  of the coarsening procedure. The default value is `0`, meaning that
  the same accuracy `eps` as for the H-Matrix assembly is used.
* ``bempp.api.global_parameters.hmat.compression_algorithm``: The compression
  algorithm that is used. The default is `aca`, which uses a modifed Adaptive
  Cross Approximation with heuristic failure detection. By choosing `dense`
  no compression is performed and all blocks are stored as dense matrices.
* ``bempp.api.global_parameters.hmat.eps``: The relative accuracy of the H-Matrix
  compression.
* ``bempp.api.global_parameters.hmat.mat_vec_parallel_levels``: The H-Matrix vector
  product is a hierarchic operation. For each node of the tree parallel tasks
  corresponding to the number of children are created. This parameter states at what
  level a serial algorithm should start to be used in order to create not too many very
  small tasks that can be inefficient (default 5).
* ``bempp.api.global_parameters.hmat.max_block_size``. For load-balancing reasons it is
  useful to restrict the maximum block size as each block is compressed using a single task
  and a very large block would occupy a single task too long. The default is 2048, but should
  be increased for very large problems.
* ``bempp.api.global_parameters.hmat.max_rank``. The maximum allowed rank for the compression of
  a block. The block compression stops at this rank even if the accuracy requirements are not
  yet achieved. The default is 30.
* ``bempp.api.global_parameters.hmat.min_block_size``. A measure for the minimum size of a block.
  The default is 20, meaning that clustering stops if fewer than 20 degrees of freedom are left
  in a cluster node.




