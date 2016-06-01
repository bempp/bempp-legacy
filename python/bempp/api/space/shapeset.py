

class Shapeset(object):
    """Access the data of the basis functions on an element."""

    def __init__(self, impl):
        self._impl = impl

    def evaluate(self, points, dof, values=True, derivatives=False):

        return self._impl.evaluate(points, dof, values, derivatives)

    @property
    def order(self):
        return self._impl.order

    @property
    def size(self):
        return self._impl.size
