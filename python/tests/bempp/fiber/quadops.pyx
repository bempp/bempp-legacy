from cython.operator import dereference as deref, preincrement
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from bempp.fiber.quadrature_options cimport QuadratureOptions
from bempp.fiber.quadrature_options cimport c_QuadratureOptions
from bempp.fiber.range_accuracy_options cimport RangeAccuracyOptions
from bempp.fiber.accuracy_options cimport AccuracyOptions, c_AccuracyOptions

ctypedef fused CanFreeze:
    QuadratureOptions
    RangeAccuracyOptions

def toggle_freeze(CanFreeze quadops):
    """ Toggles frozen status on and off """
    quadops.toggle_freeze()

def test_range_to_cpp():
    from operator import itemgetter
    should_contain = {
        0.5: (1, False), 1.5: [2, True], float('inf'): (4, False)
    }
    options = RangeAccuracyOptions(should_contain)

    cdef vector[pair[double, c_QuadratureOptions]] cvector
    options.to_cpp(cvector)

    assert cvector.size() == len(should_contain)

    cdef vector[pair[double, c_QuadratureOptions]].iterator \
            iterator = cvector.begin()
    sorted_values = sorted(should_contain.items(), key=itemgetter(0))
    for distance, values in sorted_values:
        quadop = QuadratureOptions(*values)
        if distance == float('inf'):
            assert deref(iterator).first == float('inf')
        else:
            assert abs(distance - deref(iterator).first) < 1e-8
        assert deref(iterator).second.quadratureOrder(0) == quadop(0)
        assert deref(iterator).second.quadratureOrder(1) == quadop(1)
        preincrement(iterator)

def test_accuracy_to_cpp(**kwargs):
    from itertools import product
    from numpy import arange
    acc = AccuracyOptions(**kwargs)

    cdef c_AccuracyOptions c_acc
    acc.to_cpp(c_acc)

    for distance, order in product(arange(0.05, 4, 0.1), [0, 1, 2]):
        assert c_acc.singleRegular(distance).quadratureOrder(order) \
                == acc.single_regular(distance, order)
        assert c_acc.doubleRegular(distance).quadratureOrder(order) \
                == acc.double_regular(distance, order)
    for order in [0, 1, 2]:
        assert c_acc.doubleSingular().quadratureOrder(order) \
                == acc.double_singular(order)
