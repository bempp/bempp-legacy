<%
from data_types import dtypes, scalar_cython_type
%>
from bempp.hmat.block_cluster_tree cimport c_BlockClusterTreeNode
from bempp.hmat.block_cluster_tree cimport BlockClusterTreeNode
from bempp.hmat.block_cluster_tree cimport c_BlockClusterTree
from bempp.hmat.block_cluster_tree cimport BlockClusterTree
from bempp.hmat.hmatrix_data cimport c_HMatrixData
from bempp.utils cimport shared_ptr
from bempp.utils cimport catch_exception




cdef extern from "bempp/hmat/hmatrix.hpp":
    cdef cppclass c_HMatrix "hmat::DefaultHMatrixType"[T]:
        shared_ptr[const c_BlockClusterTree] blockClusterTree() const
        shared_ptr[const c_HMatrixData[T]] data(shared_ptr[const c_BlockClusterTreeNode]&) except+catch_exception






