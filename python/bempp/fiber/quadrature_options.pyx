from collections import Sequence

cdef class QuadratureOptions:
    """ Options controlling the order of numerical quadrature

        This class can be used to specify the required order of accuracy of a
        numerical quadrature rule used to approximate a certain class of
        integrals. The order can be set either as absolute or as relative with
        respect to some default value determined by the code doing the
        integration.

        Parameters:
        -----------

        value : int
            Relative or absolute order
        is_relative : bool
            Whether the order is relative or absolute
    """
    def __init__(self, value=0, is_relative=True):
        self.__is_relative = is_relative
        if self.__is_relative:
            self.impl.setRelativeQuadratureOrder(value)
        else:
            self.impl.setAbsoluteQuadratureOrder(value)

    def __call__(self, int value=0):
        return self.impl.quadratureOrder(value)

    def __getitem__(self, key):
        return (self.value, self.is_relative)[key]
    def __len__(self):
        return 2
    def __iter__(self):
        return (self.value, self.is_relative).__iter__()

    def __str__(self):
        return str((self.value, self.is_relative))
    def __repr__(self):
        return str(self.__class__.__name__) + str(self)
    def __richcmp__(self, other, int method):
        if 2 > method or method > 3:
            raise AttributeError("No ordering operators %i" % method)
        if not hasattr(other, '__len__'):
            other = QuadratureOptions(other)
        return ((self.value, self.is_relative) == other) == (method == 2)

    property value:
        """ Relative or absolute order """
        def __get__(self):
            return self.impl.quadratureOrder(0)
        def __set__(self, value):
            if self.__is_frozen:
                raise AttributeError("Options are frozen")
            if self.__is_relative:
                self.impl.setRelativeQuadratureOrder(value)
            else:
                self.impl.setAbsoluteQuadratureOrder(value)

    property is_relative:
        """ Whether the order is relative or absolute """
        def __get__(self):
            return self.__is_relative
        def __set__(self, is_relative):
            if self.__is_frozen:
                raise AttributeError("Options are frozen")
            self.__is_relative = is_relative
            self.value = self.value # Calls setter for the property value

    cdef void toggle_freeze(self, value=None):
        self.__is_frozen = value == True if value is not None \
            else not self.__is_frozen

Sequence.register(QuadratureOptions)
