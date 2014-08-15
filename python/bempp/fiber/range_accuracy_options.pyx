from bempp.fiber.quadrature_options cimport QuadratureOptions

cdef class RangeAccuracyOptions(dict):
    """ Accuracy options with spatial range

        In practice, this is a dictionary where the keys are positive real
        numbers (defined to the accuracy
        :py:data:`RangeAccuracyOptions.__tolerance__`), and the values are
        :py:class:`QuadratureOptions`.

        All values below :py:data:`RangeAccuracyOptions.__tolerance__`, as well
        as None, `'inf'`, and float('inf')` are resolved as `float('inf')`.

        >>> from bempp.fiber import RangeAccuracyOptions
        >>> acc = RangeAccuracyOptions({0.5: (1, True), 'inf': (3, True)})
        >>> acc(0.1, 0), acc(0.6, 0)
        1, 3
        >>> acc(0.1, 2), acc(0.6, 2)
        3, 3

        In the second line above, the default order is zero and hence has no
        effect. So the first call returns the order from the range below 0.5,
        whereas the second call returns the order for the range between 0.5 and
        infinity.

        The fourth line show that effect of the default order when the accuracy
        is relative and when it is absolute (no effect).
    """
    def __cinit__(self):
        self.__tolerance__ = 1e-8
        self.__is_frozen = False

    def __init__(self, *args, **kwargs):
        from collections import Sequence
        super(RangeAccuracyOptions, self).__init__()
        # instanciate from a quadrature option
        if len(kwargs) == 0 and len(args) in [1, 2]:
            try:
                self['inf'] = args if len(args) == 2 else args[0]
                return
            except:
                pass

        self.update(dict(*args, **kwargs))

    def update(self, mapping):
        if self.__is_frozen:
            raise AttributeError("This object can no longuer be modified")
        for key, value in mapping.iteritems():
            self[key] = value

    def __setitem__(self, key, value):
        from collections import Sequence

        if self.__is_frozen:
            raise AttributeError("This object can no longuer be modified")

        cdef double c_key = self.__index(key)
        if c_key < 0:
            raise KeyError(key)
        qops = QuadratureOptions(*value) if isinstance(value, Sequence) \
                else QuadratureOptions(value)
        super(RangeAccuracyOptions, self).__setitem__(c_key, qops)

    def __getitem__(self, key):
        cdef double c_key = self.__index(key)
        if c_key < 0:
            raise KeyError(key)
        return super(RangeAccuracyOptions, self).__getitem__(c_key)

    def __delitem__(self, key):
        if self.__is_frozen:
            raise AttributeError("This object can no longuer be modified")

        cdef double c_key = self.__index(key)
        if c_key < 0:
            raise KeyError(key)
        super(RangeAccuracyOptions, self).__delitem__(c_key)

    def __contains__(self, key):
        cdef double c_key = self.__index(key)
        if c_key < 0:
            raise KeyError(key)
        return super(RangeAccuracyOptions, self).__contains__(c_key)

    cdef double __index(self, key):
        """ Transforms key to a double """
        cdef:
            double c_key
            double c_value

        try:    # Makes sure input is translatable to float
            if key is None:
                return float('inf')
            py_key = float(key)
            if py_key == float('inf'):
                return float('inf')

            c_key = py_key
            if c_key < self.__tolerance__:
                return float('inf')

            for value in super(RangeAccuracyOptions, self).iterkeys():
                c_value = value
                if c_key - self.__tolerance__ < c_value \
                        and c_key + self.__tolerance__ > c_value:
                    return c_value
            return c_key
        except: # If not, exception will be triggered outside this function
            return -1

    def fromkeys(self, *args):
        raise NotImplementedError()
    def setdefault(self, *args):
        raise NotImplementedError()

    cdef toggle_freeze(self, value=None):
        """ Allows/Disallows changing data owned by this instance """
        cdef cbool c_value = value == True if value is not None \
            else not self.__is_frozen
        self.__is_frozen = c_value
        for quadops in self.itervalues():
            (<QuadratureOptions> quadops).toggle_freeze(c_value)

    cdef void to_cpp(self,
            vector[pair[double, c_QuadratureOptions]] &quadops) except *:
        """ Converts from python to C++ values """
        cdef pair[double, c_QuadratureOptions] c_pair
        for d, quadop in self.iteritems():
            c_pair.first = <double> d
            if quadop.is_relative:
                c_pair.second.setRelativeQuadratureOrder(quadop.value)
            else:
                c_pair.second.setAbsoluteQuadratureOrder(quadop.value)
            quadops.push_back(c_pair)

    def __richcmp__(self, other, int method):
        from itertools import chain
        if 2 > method or method > 3:
            raise AttributeError("No ordering operators")

        if not isinstance(other, RangeAccuracyOptions):
            other = RangeAccuracyOptions(other)
        cdef:
            double tol = self.__tolerance__
            cbool condition = (method == 2)
        if len(self) != len(other):
            return not condition
        cdef cbool equal
        for a, b in chain(zip(self.iteritems(), other.iteritems())):
            equal =                                                     \
                (abs(a[0] - b[0]) < tol or a[0] == b[0] == float('inf')) \
                and a[1] == b[1]
            if not equal:
                return not condition
        return condition

    def __call__(self, double distance=0, int order=0):
        """ Integration accuracy

            Parameter
            ---------
            distance : double
                Distance for which to give accuracy.
            order : int
                Default order if the accuracy is relative. Ignored otherwise.
        """
        from operator import itemgetter
        for d, acc in sorted(self.iteritems(), key=itemgetter(0)):
            if distance < d + self.__tolerance__:
                return acc(order)
