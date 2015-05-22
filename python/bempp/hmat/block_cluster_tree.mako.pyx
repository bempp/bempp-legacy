from cython.operator cimport dereference as deref
from cython.operator cimport address


cdef class IndexRange:

    def __cinit__(self, size_t begin, size_t end):
        self.impl_[0] = begin
        self.impl_[1] = end

    def __init__(self, size_t begin, size_t end):
        pass

    def __dealloc__(self):
        pass

    property begin:

        def __get__(self):
            return self.impl_[0]

        def __set__(self, size_t i):
            self.impl_[0] = i

    property end:

        def __get__(self):
            return self.impl_[1]

        def __set__(self, size_t i):
            self.impl_[1] = i

cdef class BlockClusterTreeNode:

    def __cinit__(self):
        pass

    def __init__(self):
        pass


    def __dealloc__(self):
        self.impl_.reset()

    def child(self, int i):

        if i<0 or i>3:
            raise ValueError("Index must be between 0 and 3.")

        cdef BlockClusterTreeNode result = BlockClusterTreeNode()
        result.impl_.assign(deref(self.impl_).child(i))
        result._initialize()
        return result

    property children:

        def __get__(self):
            for i in range(4):
                yield self.child(i)

    property row_cluster_range:

        def __get__(self):
            cdef IndexRange result = IndexRange()

            result.impl_ = deref(deref(self.impl_).data().rowClusterTreeNode).data().indexRange
            return (result.begin,result.end)

    property column_cluster_range:

        def __get__(self):
            cdef IndexRange result = IndexRange()

            result.impl_ = deref(deref(self.impl_).data().columnClusterTreeNode).data().indexRange
            return (result.begin,result.end)

    property shape:

        def __get__(self):
            return (self.row_cluster_range[1]-self.row_cluster_range[0],
                    self.column_cluster_range[1]-self.column_cluster_range[0])

    property admissible:

        def __get__(self):
            return deref(self.impl_).data().admissible

    property is_leaf:

        def __get__(self):
            return deref(self.impl_).isLeaf()


cdef class BlockClusterTree:

    def __cinit__(self):
        pass

    def __init__(self):
        self.leaf_nodes = None

    def __dealloc__(self):
        self.impl_.reset()

    property root:

        def __get__(self):
            cdef BlockClusterTreeNode node = BlockClusterTreeNode()
            node.impl_ =  deref(self.impl_).root()
            return node

    property leaf_nodes:

        def __get__(self):

            def return_leafs(node):
                if node.is_leaf:
                    yield node
                else:
                    for child in node.children:
                        return_leafs(child)

            return_leafs(self.root)





    
