include_directories("${PYTHON_BINARY_DIR}/bempp/include")
if(WITH_TESTS)
    add_pytest(test_options.py PREFIX bempp FAKE_INIT)
endif()

add_subdirectory(grid)
add_subdirectory(fiber)
