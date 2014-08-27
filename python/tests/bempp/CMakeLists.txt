include_directories("${PYTHON_BINARY_DIR}/bempp/include")
if(WITH_TESTS)
    add_pytest(test_options.py _options.pyx
        PREFIX bempp FAKE_INIT CPP LIBRARIES libbempp)
    add_dependencies(bempp.tests._options-cython cython-headers)
endif()

add_subdirectory(grid)
add_subdirectory(fiber)
