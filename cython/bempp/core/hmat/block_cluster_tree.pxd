from bempp.core.utils.parameter_list cimport shared_ptr, c_ParameterList
from bempp.core.space.space cimport Space, c_Space
from libcpp cimport bool as cbool
from libcpp.vector cimport vector


cdef extern from "bempp/hmat/common.hpp":
    cdef cppclass c_IndexRangeType "hmat::IndexRangeType":
        c_IndexRangeType()
        size_t& operator[](size_t)

cdef extern from "bempp/hmat/bounding_box.hpp":
    cdef cppclass c_BoundingBox "hmat::BoundingBox":

        double xmin() const
        double xmax() const
        double ymin() const
        double ymax() const
        double zmin() const
        double zmax() const

        double diameter() const
        double distance(const c_BoundingBox&) const

cdef extern from "bempp/hmat/cluster_tree.hpp":
    cdef cppclass c_ClusterTreeNodeData "hmat::ClusterTreeNodeData":
        c_IndexRangeType indexRange
        c_BoundingBox boundingBox

cdef extern from "bempp/hmat/cluster_tree.hpp":
    cdef cppclass c_ClusterTreeNode "hmat::ClusterTreeNode<2>":
        const c_ClusterTreeNodeData& data() const


cdef extern from "bempp/hmat/block_cluster_tree.hpp":
    cdef cppclass c_BlockClusterTreeNodeData "hmat::BlockclusterTreeNodeData<2>":
        cbool admissible
        shared_ptr[const c_ClusterTreeNode] rowClusterTreeNode
        shared_ptr[const c_ClusterTreeNode] columnClusterTreeNode

    cdef cppclass c_BlockClusterTreeNode "hmat::BlockClusterTreeNode<2>":
        shared_ptr[const c_BlockClusterTreeNode] child(int) const
        const c_BlockClusterTreeNodeData& data() const
        cbool isLeaf() const

cdef extern from "bempp/hmat/block_cluster_tree.hpp":
    cdef cppclass c_BlockClusterTree "hmat::BlockClusterTree<2>":
        size_t rows() const
        size_t columns() const
        shared_ptr[const c_BlockClusterTreeNode] root() const

cdef extern from "bempp/assembly/hmat_interface.hpp":
    cdef shared_ptr[const c_BlockClusterTree] c_generateBlockClusterTree "Bempp::generateBlockClusterTree" [BASIS](
            const c_Space[BASIS]&,
            const c_Space[BASIS]&,
            const c_ParameterList&)

cdef class BlockClusterTreeNode:
    cdef shared_ptr[const c_BlockClusterTreeNode] impl_

cdef class BlockClusterTree:
    cdef shared_ptr[const c_BlockClusterTree] impl_

cdef class BoundingBox:
    cdef c_BoundingBox impl_





