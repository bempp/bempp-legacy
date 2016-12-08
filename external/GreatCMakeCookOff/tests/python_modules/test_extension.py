def test_imports():
    """ Extension contains expected attributes """
    import extension
    assert 'Error' in dir(extension)
    assert 'error_out' in dir(extension)
    assert 'meaning_of_life' in dir(extension)


def test_Error():
    """ Error object derives from Exception """
    from extension import Error
    assert Error.__base__ is Exception


def test_error_out():
    """ Creates error via call to C """
    from pytest import raises
    from extension import error_out
    with raises(Exception):
        error_out("Hello, world")


def test_meaning_of_life():
    """ Function in other file is included in extension """
    from extension import meaning_of_life
    assert meaning_of_life() == 42
