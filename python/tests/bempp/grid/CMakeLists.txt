if(WITH_TESTS)
    add_pytest(test_grid.py test_grid_factory.py 
               PREFIX bempp.grid FAKE_INIT)
endif()
