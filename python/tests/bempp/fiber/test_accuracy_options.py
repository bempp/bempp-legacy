from py.test import mark
from bempp.fiber import AccuracyOptions

inputs = [
    {
        'single_regular': (1, False),
        'double_regular': (2, True),
        'double_singular': (3, True)
    },
    {
        'single_regular': 1,
        'double_regular': {1: (2, True), 'inf': 5},
        'double_singular': 2,
    },
    {
        'single_regular': {0.1: (2, True), 0.3: 5, 'inf': 5},
        'double_regular': {1: (2, True), 0.2: (4, False), 'inf': 5},
        'double_singular': 2,
    }
]


@mark.parametrize("kwargs", inputs)
def test_instantiation(kwargs):
    acc = AccuracyOptions(**kwargs)
    assert acc.double_singular == kwargs['double_singular']
    assert acc.single_regular == kwargs['single_regular']
    assert acc.double_regular == kwargs['double_regular']


@mark.parametrize("kwargs", inputs)
def test_equality(kwargs):
    acc = AccuracyOptions(**kwargs)

    for equal in (acc, AccuracyOptions(**kwargs)):
        assert acc == equal
        assert not (acc != equal)


@mark.parametrize("kwargs", inputs)
def test_to_cpp(kwargs):
    from bempp.fiber.tests.quadops import test_accuracy_to_cpp
    test_accuracy_to_cpp(**kwargs)
