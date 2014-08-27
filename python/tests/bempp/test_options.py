""" Test general BEM++ option class

    For simplicity, only test a subset of options one per type more or less.
"""
from bempp.options import Options
from py.test import mark, fixture


@fixture
def from_kwargs():
    return Options(
        eps=1e-8, maximum_block_size=100, recompress=False,
        aca_assembly_mode="local"
    )


@fixture
def from_set():
    options = Options()

    options.eps = 1e-8
    options.maximum_block_size = 100
    options.recompress = False,
    options.aca_assembly_mode = "local"

    return options


@mark.parametrize("instantiater", [from_kwargs, from_set])
def test_options(instantiater):
    options = instantiater()

    assert abs(options.eta - 1e-4) < 1e-8
    assert abs(options.eps - 1e-8) < 1e-8
    assert options.maximum_block_size == 100
    assert options.minimum_block_size == 16
    assert options.maximum_rank in [2**32-2, 2**64-2]
    assert options.reaction_to_unsupported_mode == "warning"
    assert options.aca_assembly_mode == "local"
