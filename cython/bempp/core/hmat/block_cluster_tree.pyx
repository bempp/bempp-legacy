# cython: embedsignature=True

from cython.operator cimport dereference as deref
from cython.operator cimport address

from bempp.core.space.space cimport Space, c_Space
from bempp.core.utils.parameter_list cimport ParameterList

cdef class BlockClusterTreeNode:
    """Interface to a single node of the block cluster tree."""

    def __cinit__(self):
        pass

    def __init__(self):
        pass


    def __dealloc__(self):
        self.impl_.reset()

    def child(self, int i):
        """Access to a child node by index 0 <= i < 4."""

        if i<0 or i>3:
            raise ValueError("Index must be between 0 and 3.")

        cdef BlockClusterTreeNode result = BlockClusterTreeNode()
        result.impl_.assign(deref(self.impl_).child(i))
        return result

    property children:
        """Return an iterator over all child nodes."""

        def __get__(self):
            for i in range(4):
                yield self.child(i)

    property row_cluster_range:
        """Return the range of the associated row cluster node."""

        def __get__(self):
            s = deref(deref(self.impl_).data().rowClusterTreeNode).data().indexRange[0]
            e = deref(deref(self.impl_).data().rowClusterTreeNode).data().indexRange[1]
            return (s, e)

    property column_cluster_range:
        """Return the range of the associated column cluster node."""

        def __get__(self):
            s = deref(deref(self.impl_).data().columnClusterTreeNode).data().indexRange[0]
            e = deref(deref(self.impl_).data().columnClusterTreeNode).data().indexRange[1]
            return (s, e)

    property shape:
        """Return the shape of the block represented by this node."""

        def __get__(self):
            return (self.row_cluster_range[1]-self.row_cluster_range[0],
                    self.column_cluster_range[1]-self.column_cluster_range[0])

    property admissible:
        """Return True if the node is admissible, otherwise False."""

        def __get__(self):
            return deref(self.impl_).data().admissible

    property is_leaf:
        """Return True if the node is a leaf, otherwise False."""

        def __get__(self):
            return deref(self.impl_).isLeaf()


cdef class BlockClusterTree:
    """Interface to the block cluster tree structre of an H-Matrix."""

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def plot(self,file_name = 'block_cluster_tree.png', display = True, width = 2000, height= 2000,
            delete=True):
        """Plot a given block cluster tree into a file. If `display` is True the plot is also displayed.

        Parameters
        ----------
        file_name : string
            name of the output file.
        display : bool
            If True also display the plot.
        width : int
            Width in pixels.
        height : int
            Height in pixels.
        delete : bool
            If True image is deleted after displaying it.

        Notes
        -----
        This method requires that PyQt4 is installed.

        """

        from PyQt4 import QtGui

        rows = self.root.shape[0]
        cols = self.root.shape[1]

        row_scale = (1.*height)/rows
        col_scale = (1.*width)/cols

        img = QtGui.QImage(width, height,QtGui.QImage.Format_RGB32)
        qp = QtGui.QPainter(img)
        color = QtGui.QColor()

        qp.setPen(QtGui.QColor(0,0,0))
        for node in self.leaf_nodes:
            if node.admissible:
                qp.setBrush(QtGui.QColor(0,255,0))
            else:
                qp.setBrush(QtGui.QColor(255,0,0))
            qp.drawRect(node.column_cluster_range[0]*col_scale,
                        node.row_cluster_range[0]*row_scale,
                        node.shape[1]*col_scale,
                        node.shape[0]*row_scale)

        qp.end()
        img.save(file_name)
        
        if display:
            from matplotlib import image as mpimage
            from matplotlib import pyplot as plt

            img = mpimage.imread(file_name)
            plt.imshow(img,extent=[0,cols,rows,0])
            plt.show()
            if delete:
                import os
                os.remove(file_name)

    property root:
        """Return the root node of the tree."""

        def __get__(self):
            cdef BlockClusterTreeNode node = BlockClusterTreeNode()
            node.impl_ =  deref(self.impl_).root()
            return node

    property leaf_nodes:
        """Return an iterator over the leaf nodes."""

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
        ParameterList parameter_list = None):
    """Create a block cluster tree for a given pair of spaces."""

    from bempp.api import global_parameters

    cdef BlockClusterTree result = BlockClusterTree()
    cdef ParameterList parameters

    if parameter_list is not None:
        parameters = parameter_list
    else:
        parameters = global_parameters

    result.impl_.assign(c_generateBlockClusterTree[double](
        deref(test_space.impl_),
        deref(trial_space.impl_),
        deref(parameters.impl_)))
    return result




    
