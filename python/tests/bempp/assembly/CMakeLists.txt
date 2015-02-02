if (WITH_TESTS)
    add_pytest(test_grid_function.py PREFIX bempp.assembly FAKE_INIT)
    add_pytest(test_discrete_boundary_operator.py PREFIX bempp.assembly FAKE_INIT)
endif()

