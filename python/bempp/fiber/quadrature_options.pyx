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

    def order(self, value=0):
        return self.impl.quadratureOrder(value)

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
