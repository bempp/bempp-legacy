def test_imports():
    from extension import structure
    assert 'Structure' in dir(structure)


def test_meaning_of_life():
    from py.test import raises
    from extension.structure import Structure
    structure = Structure()
    assert structure.meaning_of_life == 42
    structure.meaning_of_life = 33
    assert structure.meaning_of_life == 33
    structure.meaning_of_life = 34.3
    assert structure.meaning_of_life == 34
    structure.meaning_of_life = "55"
    assert structure.meaning_of_life == 55
    with raises(ValueError):
        structure.meaning_of_life = "hello"


def test_message():
    from sys import version
    from py.test import raises
    from extension.structure import Structure
    structure = Structure()
    if int(version[0]) < 3:
        assert structure.message == "And everything"
    else:
        assert structure.message == bytes("And everything", 'utf-8')
    with raises(AttributeError):
        structure.message = "hello"
