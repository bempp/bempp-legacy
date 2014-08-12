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
    assert ops.order(0) == 2
    assert ops.order(1) == 3

    # Change value
    ops.value = 55
    assert ops.value == 55
    assert ops.is_relative
    assert ops.order(-1) == 54

    # Change to absolute
    ops.is_relative = False
    assert ops.value == 55
    assert not ops.is_relative
    assert ops.order(-1) == 55


def test_absolute():
    from bempp.fiber import QuadratureOptions
    ops = QuadratureOptions(2, False)
    assert ops.value == 2
    assert not ops.is_relative
    assert ops.order(0) == 2
    assert ops.order(-1) == 2

    # Change value
    ops.value = 55
    assert ops.value == 55
    assert not ops.is_relative
    assert ops.order(1) == 55

    # Change to relative
    ops.is_relative = True
    assert ops.value == 55
    assert ops.is_relative
    assert ops.order(-1) == 54
