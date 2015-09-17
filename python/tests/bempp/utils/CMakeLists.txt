if(WITH_TESTS)
    add_pytest(test_blocked_linear_operator.py 
               test_eigen_conversion.py PREFIX bempp.utils FAKE_INIT)
    add_pytest(test_pickle_global_options.py PREFIX bempp.utils FAKE_INIT)
endif()
