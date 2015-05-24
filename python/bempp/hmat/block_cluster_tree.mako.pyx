from cython.operator cimport dereference as deref
from cython.operator cimport address

from bempp.space.space cimport Space, c_Space, _py_get_space_ptr

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
        return result

    property children:

        def __get__(self):
            for i in range(4):
                yield self.child(i)

    property row_cluster_range:

        def __get__(self):
            cdef IndexRange result = IndexRange(0,0)

            result.impl_ = deref(deref(self.impl_).data().rowClusterTreeNode).data().indexRange
            return (result.begin,result.end)

    property column_cluster_range:

        def __get__(self):
            cdef IndexRange result = IndexRange(0,0)

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
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def plot(self):

        from matplotlib import pyplot as plt
        from matplotlib.patches import Rectangle
        from matplotlib.collections import PatchCollection

        rows = self.root.shape[0]
        cols = self.root.shape[1]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_xlim([0,cols])
        ax.set_ylim([0,rows])
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.xaxis.tick_top()
        ax.yaxis.tick_left()

        for node in self.leaf_nodes:
            if node.admissible:
                color = 'green'
            else:
                color = 'red'
            ax.add_patch(Rectangle(
                    (node.column_cluster_range[0],
                     node.row_cluster_range[0]),
                     node.shape[1],
                     node.shape[0],
                     color=color))

        plt.show()




    property root:

        def __get__(self):
            cdef BlockClusterTreeNode node = BlockClusterTreeNode()
            node.impl_ =  deref(self.impl_).root()
            return node

    property leaf_nodes:

        def __get__(self):

            leafs = []

            def iterate_leafs(node):
                if node.is_leaf:
                    leafs.append(node)
                else:
                    for child in node.children:
                        iterate_leafs(child)

            iterate_leafs(self.root)

            for node in leafs:
                yield node


def generate_block_cluster_tree(Space test_space, Space trial_space, 
        int min_block_size, 
        int max_block_size, 
        double eta):

    cdef BlockClusterTree result = BlockClusterTree()
    result.impl_.assign(c_generateBlockClusterTree[double](
        deref(_py_get_space_ptr[double](test_space.impl_)),
        deref(_py_get_space_ptr[double](trial_space.impl_)),
        min_block_size,
        max_block_size,
        eta))
    return result




    
