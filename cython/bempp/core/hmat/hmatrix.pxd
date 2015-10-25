from bempp.core.hmat.block_cluster_tree cimport c_BlockClusterTreeNode
from bempp.core.hmat.block_cluster_tree cimport BlockClusterTreeNode
from bempp.core.hmat.block_cluster_tree cimport c_BlockClusterTree
from bempp.core.hmat.block_cluster_tree cimport BlockClusterTree
from bempp.core.hmat.hmatrix_data cimport c_HMatrixData
from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport catch_exception

cdef extern from "bempp/hmat/hmatrix.hpp":
    cdef cppclass c_HMatrix "hmat::DefaultHMatrixType"[T]:
        shared_ptr[const c_BlockClusterTree] blockClusterTree() const
        shared_ptr[const c_HMatrixData[T]] data(shared_ptr[const c_BlockClusterTreeNode]&) except+catch_exception
        double frobeniusNorm() const
        int numberOfDenseBlocks() const
        int numberOfLowRankBlocks() const
        int numberOfBlocks() const
        double memSizeKb();






