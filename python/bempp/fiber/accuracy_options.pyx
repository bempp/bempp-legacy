cdef class AccuracyOptions:
    """ Quadrature accuracy options

        Parameters:
        -----------
        single_regular : RangeAccuracyOptions
            Accuracy of the integration of regular functions on single elements
        double_regular : RangeAccuracyOptions
            Accuracy of the integration of regular functions on pair of
            elements
        double_singular : QuadratureOptions
            Accuracy of the integration of singular functions on pair of
            elements
    """
    def __cinit__(self, single_regular, double_regular, double_singular):
        self.single_regular = RangeAccuracyOptions(single_regular)
        self.double_regular = RangeAccuracyOptions(double_regular)
        self.double_singular = double_singular
    cdef toggle_freeze(self, value=None):
        self.single_regular.toggle_freeze(value)
        self.double_regular.toggle_freeze(value)
        self.double_singular.toggle_freeze(value)

    property single_regular:
        """ Integration accuracy of regular functions on single elements """
        def __get__(self):
            return self.__single_regular
        def __set__(self, value):
            self.__single_regular = RangeAccuracyOptions(value)

    property double_regular:
        """ Integration accuracy of regular functions on pair of elements """
        def __get__(self):
            return self.__double_regular
        def __set__(self, value):
            self.__double_regular = RangeAccuracyOptions(value)

    property double_singular:
        """ Integration accuracy of singular functions on pair of elements """
        def __get__(self):
            return self.__double_singular
        def __set__(self, value):
            from collections import Sequence
            if isinstance(value, Sequence):
                self.__double_singular = QuadratureOptions(*value)
            else:
                self.__double_singular = QuadratureOptions(value)

    def __str__(self):
        return "{0.__class__.__name__}:\n"                \
                "   single regular: {0.single_regular}\n" \
                "   double regular: {0.double_regular}\n" \
                "   double singular: {0.double_singular}".format(self)
    def __repr__(self):
        return "{0.__class__.__name__}("                \
                "single_regular={0.single_regular!r}, " \
                "double_regular={0.double_regular!r}, " \
                "double_singular={0.double_singular!r})".format(self)

    def __richcmp__(self, other, int method):
        from collections import Mapping, Sequence
        from itertools import chain
        if 2 > method or method > 3:
            raise AttributeError("No ordering operators")

        if isinstance(other, Mapping):
            other = other['single_regular'], \
                other['double_regular'], \
                other['double_singular']
        elif isinstance(other, AccuracyOptions):
            other = other.single_regular, \
                other.double_regular, \
                other.double_singular

        this = self.single_regular, self.double_regular, self.double_singular
        return (this == other) == (method == 2)

    cdef void to_cpp(self, c_AccuracyOptions &c_options) except*:

        c_options.setDoubleSingular(self.__double_singular.value,
            self.__double_singular.is_relative)

        cdef vector[pair[double, c_QuadratureOptions]] c_regular
        self.__single_regular.to_cpp(c_regular)
        c_options.setSingleRegular(c_regular)

        c_regular.clear()
        self.__double_regular.to_cpp(c_regular)
        c_options.setDoubleRegular(c_regular)
