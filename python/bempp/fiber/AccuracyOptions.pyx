from bempp.fiber.quadrature_options cimport QuadratureOptions

def __convert_to_options(options, must_be_set = False):
    if isinstance(options, QuadratureOptions):
        return set(options) if must_be_set else options
    if len(options) == 2:
        try:
            result = QuadratureOptions(*options)
        except: pass
        else:
            return set(options) if must_be_set else options
    for 



cdef class AccuracyOptions:
    def __cinit__(self,
            all_or_single_regular,
            double_regular=None,
            double_singular=None):
        if double_regular is None and double_singular is None:
            self.double_singular = all_or_single_regular
        elif double_regular is not None and double_singular is not None:
            self.single_regular = all_or_single_regular
            self.double_regular = double_regular
            self.double_singular = double_singular
        else:
            msg = "Input must either specify accuracy for one or all"
            raise ArgumentError(msg)
