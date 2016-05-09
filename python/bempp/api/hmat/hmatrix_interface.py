"""Interface functions to access HMatrix properties."""


def block_cluster_tree(discrete_operator):
    """Return the block cluster tree for a HMatrix operator."""
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import block_cluster_tree_ext
    return block_cluster_tree_ext(discrete_operator._impl)


def number_of_dense_blocks(discrete_operator):
    """Return the number of dense blocks for a HMatrix operator."""
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import number_of_dense_blocks_ext
    return number_of_dense_blocks_ext(discrete_operator._impl)


def number_of_low_rank_blocks(discrete_operator):
    """Return the number of low rank blocks for a HMatrix operator."""
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import number_of_low_rank_blocks_ext
    return number_of_low_rank_blocks_ext(discrete_operator._impl)


def number_of_blocks(discrete_operator):
    """Return the number of blocks for a HMatrix operator."""
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import number_of_blocks_ext
    return number_of_blocks_ext(discrete_operator._impl)


def mem_size(discrete_operator):
    """Return the memory size in kb for a HMatrix operator."""
    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import mem_size_ext
    return mem_size_ext(discrete_operator._impl)


def data_block(discrete_operator, block_cluster_tree_node):
    """Return the data block associated with a block cluster tree node."""

    from bempp.api.assembly.discrete_boundary_operator import \
        GeneralNonlocalDiscreteBoundaryOperator

    if not isinstance(discrete_operator, GeneralNonlocalDiscreteBoundaryOperator):
        raise ValueError("discrete operator is not an HMatrix operator.")
    from bempp.core.hmat.hmatrix_interface import data_block_ext
    return data_block_ext(discrete_operator._impl, block_cluster_tree_node)


def compression_rate(discrete_operator):
    """Return the compression rate for a HMatrix operator."""

    mem = mem_size(discrete_operator)
    total = 8 if discrete_operator.dtype == 'float64' else 16
    total *= discrete_operator.shape[0] * \
        discrete_operator.shape[1] / (1.0 * 1024)
    return mem / total
