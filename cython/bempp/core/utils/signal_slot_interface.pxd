
cdef extern from "<boost/signals2/signal.hpp>" namespace "boost":
    cdef cppclass Connection "boost::signals2::connection":
        disconnect()

cdef extern from "bempp/python/slot_interface.hpp" namespace "Bempp":
    cdef cppclass SlotInterface:
        SlotInterface(object)

