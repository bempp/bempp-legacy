if (WITH_TESTS)
    add_pytest(test_p1_coupling.py PREFIX bempp.fenics_interface FAKE_INIT)
endif()

