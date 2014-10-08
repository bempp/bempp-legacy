def test_default_creation():
    from bempp.fiber import QuadratureOptions
    ops = QuadratureOptions()
    assert ops.value == 0
    assert ops.is_relative


def test_relative():
    from bempp.fiber import QuadratureOptions
    ops = QuadratureOptions(2, True)
    assert ops.value == 2
    assert ops.is_relative
    assert ops(0) == 2
    assert ops(1) == 3

    # Change value
    ops.value = 55
    assert ops.value == 55
    assert ops.is_relative
    assert ops(-1) == 54

    # Change to absolute
    ops.is_relative = False
    assert ops.value == 55
    assert not ops.is_relative
    assert ops(-1) == 55


def test_absolute():
    from bempp.fiber import QuadratureOptions
    ops = QuadratureOptions(2, False)
    assert ops.value == 2
    assert not ops.is_relative
    assert ops(0) == 2
    assert ops(-1) == 2

    # Change value
    ops.value = 55
    assert ops.value == 55
    assert not ops.is_relative
    assert ops(1) == 55

    # Change to relative
    ops.is_relative = True
    assert ops.value == 55
    assert ops.is_relative
    assert ops(-1) == 54


def test_frozen():
    from py.test import raises
    from bempp.fiber import QuadratureOptions
    from bempp.fiber.tests.quadops import toggle_freeze
    ops = QuadratureOptions(2, False)
    assert ops.value == 2
    assert not ops.is_relative

    ops.value = 6
    ops.is_relative = True
    assert ops.value == 6
    assert ops.is_relative

    toggle_freeze(ops)
    with raises(AttributeError):
        ops.value = 5
    with raises(AttributeError):
        ops.is_relative = True


def test_equality():
    from bempp.fiber import QuadratureOptions

    quadops = QuadratureOptions(1, True)
    is_equal = [quadops, QuadratureOptions(1, True), (1, True), 1]
    is_not_equal = [
        QuadratureOptions(2, True), (2, True),
        QuadratureOptions(1, False), (1, False),
        QuadratureOptions(2, False), (2, False),
        3
    ]

    for b in is_equal:
        assert quadops == b
        assert not (quadops != b)
    for b in is_not_equal:
        assert not (quadops == b)
        assert quadops != b


def test_ordering_operators_not_defined():
    from py.test import raises
    from bempp.fiber import QuadratureOptions

    quadops = QuadratureOptions(1, False)
    with raises(AttributeError):
        quadops < (1, False)
    with raises(AttributeError):
        quadops <= (1, False)
    with raises(AttributeError):
        quadops >= (1, False)
    with raises(AttributeError):
        quadops > (1, False)


def test_pass_as_varargs():
    from bempp.fiber import QuadratureOptions
    expected = (1, False)
    quadops = QuadratureOptions(*expected)

    def tester(*args):
        assert args[0] == expected[0]
        assert args[1] == expected[1]

    tester(*quadops)


def test_QuadratureOptions_is_sequence():
    from collections import Sequence
    from bempp.fiber import QuadratureOptions

    quadops = QuadratureOptions()
    assert isinstance(quadops, Sequence)
