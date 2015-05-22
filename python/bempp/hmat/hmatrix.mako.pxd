<%
from data_types import dtypes, scalar_cython_type
%>
from bempp.hmat.block_cluster_tree cimport c_BlockClusterTreeNode
from bempp.hmat.block_cluster_tree cimport BlockClusterTreeNode
from bempp.hmat.block_cluster_tree cimport c_BlockClusterTree
from bempp.hmat.block_cluster_tree cimport BlockClusterTree
from bempp.utils cimport shared_ptr




cdef extern from "bempp/hmat/hmatrix.hpp":
    cdef cppclass c_HMatrix "hmat::DefaultHMatrixType"[T]:
        shared_ptr[const c_BlockClusterTree] blockClusterTree() const





