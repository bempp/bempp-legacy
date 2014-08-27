from bempp.options cimport Options, AcaOptions, \
        WARNING, LOCAL_ASSEMBLY
def test_cconversion(Options options):
    cdef AcaOptions c_options
    options.to_aca_options(&c_options)

    assert abs(c_options.eta - 1e-4) < 1e-8
    assert abs(c_options.eps - 1e-8) < 1e-8
    assert c_options.maximumBlockSize == 100
    assert c_options.minimumBlockSize == 16
    assert c_options.maximumRank in [2**32-2, 2**64-2]
    assert c_options.reactionToUnsupportedMode == WARNING
    assert c_options.mode == LOCAL_ASSEMBLY
