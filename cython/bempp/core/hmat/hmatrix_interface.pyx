from bempp.core.hmat.block_cluster_tree cimport c_BlockClusterTreeNode
from bempp.core.hmat.block_cluster_tree cimport BlockClusterTreeNode
from bempp.core.hmat.block_cluster_tree cimport c_BlockClusterTree
from bempp.core.hmat.block_cluster_tree cimport BlockClusterTree
from bempp.core.hmat.hmatrix_data cimport c_HMatrixData
from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport catch_exception
from bempp.core.utils cimport complex_double
from bempp.core.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.core.assembly.discrete_boundary_operator cimport RealDiscreteBoundaryOperator
from bempp.core.assembly.discrete_boundary_operator cimport ComplexDiscreteBoundaryOperator
from cython.operator cimport dereference as deref


cdef extern from "bempp/hmat/hmatrix.hpp":
    cdef cppclass c_HMatrix "hmat::DefaultHMatrixType"[T]:
        shared_ptr[const c_BlockClusterTree] blockClusterTree() const
        shared_ptr[const c_HMatrixData[T]] data(shared_ptr[const c_BlockClusterTreeNode]&) except+catch_exception
        double frobeniusNorm() const
        int numberOfDenseBlocks() const
        int numberOfLowRankBlocks() const
        int numberOfBlocks() const
        double memSizeKb();

cdef extern from "bempp/assembly/discrete_hmat_boundary_operator.hpp" namespace "Bempp":
    cdef shared_ptr[const c_HMatrix[T]] castToHMatrix[T](
            const shared_ptr[const c_DiscreteBoundaryOperator[T]]&) except+catch_exception


def block_cluster_tree_ext(discrete_operator):
    """Return the block cluster tree of a discrete operator."""

    cdef BlockClusterTree block_cluster_tree = BlockClusterTree()
    try:
        if discrete_operator.dtype == 'float64':
            block_cluster_tree.impl_.assign(
                    deref(castToHMatrix[double]((<RealDiscreteBoundaryOperator>discrete_operator).impl_)).blockClusterTree())
        else:
            block_cluster_tree.impl_.assign(
                    deref(castToHMatrix[complex_double]((<ComplexDiscreteBoundaryOperator>discrete_operator).impl_)).blockClusterTree())
    except:
        raise ValueError("discrete_operator does not seem to be a valid HMatrix.")
    return block_cluster_tree


def number_of_dense_blocks_ext(discrete_operator):
    """Return the number of dense blocks of a discrete operator."""

    try:
        if discrete_operator.dtype == 'float64':
            return deref(castToHMatrix[double]((
                <RealDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfDenseBlocks()
        else:
            return deref(castToHMatrix[complex_double]((
                <ComplexDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfDenseBlocks()
    except:
        raise ValueError("discrete_operator does not seem to be a valid HMatrix.")

def number_of_low_rank_blocks_ext(discrete_operator):
    """Return the number of low rank blocks of a discrete operator."""

    try:
        if discrete_operator.dtype == 'float64':
            return deref(castToHMatrix[double]((
                <RealDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfLowRankBlocks()
        else:
            return deref(castToHMatrix[complex_double]((
                <ComplexDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfLowRankBlocks()
    except:
        raise ValueError("discrete_operator does not seem to be a valid HMatrix.")

def number_of_blocks_ext(discrete_operator):
    """Return the number of blocks of a discrete operator."""

    try:
        if discrete_operator.dtype == 'float64':
            return deref(castToHMatrix[double]((
                <RealDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfBlocks()
        else:
            return deref(castToHMatrix[complex_double]((
                <ComplexDiscreteBoundaryOperator>discrete_operator).impl_)).numberOfBlocks()
    except:
        raise ValueError("discrete_operator does not seem to be a valid HMatrix.")

def mem_size_ext(discrete_operator):
    """Return the memory size in kb."""

    try:
        if discrete_operator.dtype == 'float64':
            return deref(castToHMatrix[double]((
                <RealDiscreteBoundaryOperator>discrete_operator).impl_)).memSizeKb()
        else:
            return deref(castToHMatrix[complex_double]((
                <ComplexDiscreteBoundaryOperator>discrete_operator).impl_)).memSizeKb()
    except:
        raise ValueError("discrete_operator does not seem to be a valid HMatrix.")


