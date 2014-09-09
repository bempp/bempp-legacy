from py.test import fixture
from bempp.nested_attributes import NestedAttributes


@fixture
def nested():
    return {
        ('this', ): 'first',
        ('that', ): 'second',
        ('there', 'a'): 'third',
        ('there', 'b'): 'fourth'
    }


def test_dir(nested):
    attr = NestedAttributes(nested)
    actual = dir(attr)
    expected = [k[0] for k in nested]
    assert set(actual) == set(expected)


def test_nested_dir(nested):
    attr = NestedAttributes(nested)
    actual = dir(attr.there)
    expected = [k[1] for k in nested if len(k) > 1]
    assert set(actual) == set(expected)


def test_access(nested):
    attr = NestedAttributes(nested)

    notnested = {k[0]: v for k, v in nested.iteritems() if len(k) == 1}
    for k, v in notnested.iteritems():
        assert getattr(attr, k) == v


def test_nested_access(nested):
    attr = NestedAttributes(nested).there

    notnested = {k[1]: v for k, v in nested.iteritems() if len(k) == 2}
    for k, v in notnested.iteritems():
        assert getattr(attr, k) == v


def test_initialization_failures():
    from py.test import raises

    with raises(TypeError):
        NestedAttributes({'this': 1})

    with raises(ValueError):
        NestedAttributes({('this', ): 1, ('this', 'that'): 2})
