file(GLOB_RECURSE PYTHON_TEST_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
# Add python unit tests into ctest framework
foreach(source ${PYTHON_TEST_SOURCES})
    if(NOT "${source}" STREQUAL "run_tests.py")
        get_filename_component(testname "${source}" NAME_WE)
        set(source "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
        add_test(NAME ${testname}
            COMMAND ${LOCAL_PYTHON_EXECUTABLE} -m pytest ${source}
        )
        set_tests_properties(${testname} PROPERTIES LABELS "python;unit")
    endif()
endforeach()

if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    include(TargetCopyFiles)
    add_custom_target(copy_python_tests ALL)
    add_copy_files(copy_python_tests
        FILES ${PYTHON_TEST_SOURCES}
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )
endif ()

# Create a fake init file to make tests subdirectory a package
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/__init__.py")
install_python(FILES
    ${PYTHON_TEST_SOURCES}
    "${CMAKE_CURRENT_BINARY_DIR}/__init__.py"
    DESTINATION bempp/tests
)
