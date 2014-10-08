if(WITH_TESTS)
    include_directories(${PROJECT_SOURCE_DIR}/python)
    add_pytest(
        test_quadrature_options.py
        test_range_accuracy_options.py
        test_accuracy_options.py
        quadops.pyx
        CPP
        PREFIX bempp.fiber
        FAKE_INIT
        LIBRARIES libbempp
    )
endif()
