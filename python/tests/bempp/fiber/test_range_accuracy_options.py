from py.test import mark
from bempp.fiber import RangeAccuracyOptions, QuadratureOptions


def test_contains():
    """ Check indexing is in interval +/- tolerance """
    options = RangeAccuracyOptions()
    assert 0.5 not in options
    assert '1' not in options

    options[0.5] = (1, True)
    assert 0.5 in options
    assert 0.5 + 0.5 * options.__tolerance__ in options
    assert 0.5 - 0.5 * options.__tolerance__ in options
    assert 0.5 + 1.5 * options.__tolerance__ not in options
    assert 0.5 - 1.5 * options.__tolerance__ not in options


def test_infinity_index():
    """ Check indexing with infinity works """
    # Create empty options
    options = RangeAccuracyOptions()

    # Should have only on option, 'inf'. Tested more fully elsewhere.
    assert len(options) == 1 and 'inf' in options
    options[None] = (1, True)
    assert len(options) == 1 and 'inf' in options
    options['inf'] = (1, True)
    assert len(options) == 1 and 'inf' in options
    options[-1] = (1, True)
    assert len(options) == 1 and 'inf' in options


@mark.parametrize("kwargs", [
    {0.5: (2, True), 0.1: (1, False)},
    [(0.5, (2, True)), (0.1, (1, False))],
    {0.5: QuadratureOptions(2, True)},
    {0.5: QuadratureOptions(2, True), 0.1: (1, False)},
    {'inf': (1, False)},
    {None: (1, False)},
    {0.5: 1}
])
def test_instantiation(kwargs):
    options = RangeAccuracyOptions(kwargs)
    iterator = kwargs.iteritems() if hasattr(kwargs, 'iteritems') else kwargs
    for key, value in iterator:
        assert key in options
        assert isinstance(options[key], QuadratureOptions)
        assert options[key] == value
    assert len(options) == len(kwargs)


def test_equality():
    kwargs = {0.5: (2, True), 0.1: 4, 'inf': (1, False)}
    options = RangeAccuracyOptions(kwargs)

    equals = [options, kwargs, RangeAccuracyOptions(kwargs)]
    nequals = [
        {},
        {0.5: (2, True)},
        {0.6: (2, True), 0.1: 4, 'inf': (1, False)},
        {0.5: (3, True), 0.1: 4, 'inf': (1, False)},
        {0.5: (2, False), 0.1: 4, 'inf': (1, False)},
        {0.5: (2, True), 0.1: 4, 0.6: 3, 'inf': (1, False)}
    ]
    for equal in equals:
        assert options == equal
        assert not (options != equal)
    for nequal in nequals:
        assert not (options == nequal)
        assert options != nequal


@mark.parametrize("args", [
    QuadratureOptions(1, False),
    (1, False),
    1
])
def test_instantiation_from_QuadratureOptions(args):
    options = RangeAccuracyOptions(args)
    assert len(options) == 1
    assert 'inf' in options
    assert options['inf'] == args


def test_getsetitems():
    options = RangeAccuracyOptions({0.5: (1, False), 1.5: (2, True)})
    # Check state after creation
    assert 0.5 in options and 1.5 in options
    assert options[0.5] == (1, False) and options[1.5] == (2, True)
    assert 2 not in options and 'inf' not in options

    options[0.5].value = 4
    options[1.5] = 3, False
    options[2] = 5
    options['inf'] = 6, False

    should_contain = {
        0.5: (4, False),
        1.5: (3, False),
        2: (5, True),
        'inf': (6, False),
    }
    for key, value in should_contain.iteritems():
        assert key in options
        assert value == options[key]


def test_unimplemented_methods():
    """ Methods that don't really make sense here """
    from py.test import raises

    options = RangeAccuracyOptions()
    with raises(NotImplementedError):
        options.fromkeys('1')
    with raises(NotImplementedError):
        options.setdefault('1')


def test_frozen():
    from py.test import raises
    from bempp.fiber.tests.quadops import toggle_freeze

    options = RangeAccuracyOptions({0.5: (1, False), 1.5: (2, True)})
    assert options[0.5].value == 1
    options[0.5].value = 2
    assert options[0.5].value == 2
    options[0.3] = 6
    assert options[0.3].value == 6

    print repr(options), type(options)
    toggle_freeze(options)
    with raises(AttributeError):
        options[0.2] = 1, False
    with raises(AttributeError):
        options[0.3].value = 5


def test_failing_index():
    """ Checks that exceptions are thrown

        Main indexing function is pure C and does not throw itself.
        So, we want to make sure that exceptions are still thrown outside of
        it.
    """
    from py.test import raises

    options = RangeAccuracyOptions({0.5: (1, False), 1.5: (2, True)})
    with raises(KeyError):
        options['a']
    with raises(KeyError):
        options[1j]


def test_cpp_conversion():
    from bempp.fiber.tests.quadops import test_range_to_cpp
    test_range_to_cpp()


def test_order():
    acc = RangeAccuracyOptions({0.1: (1, 0), 0.5: (2, 1), 'inf': (3, 1)})
    inputs = {
        # distance: { default order: result from calling acc  }
        0: {0: 1, 1: 1},
        0.05: {0: 1, 1: 1},
        0.25: {0: 2, 1: 3},
        0.55: {0: 3, 1: 4},
    }
    for distance, orders in inputs.iteritems():
        for order, expected in orders.iteritems():
            assert acc(distance, order) == expected


def test_default_instantiation():
    acc = RangeAccuracyOptions()
    assert len(acc) == 1
    assert 'inf' in acc
    assert acc(0) == 0
    assert acc(1) == 0
