file(GLOB_RECURSE PYTHON_TEST_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
if (NOT (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR}))
    foreach (py_file ${PYTHON_TEST_SOURCES})
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
            COMMAND ${CMAKE_COMMAND} copy -E
                ${py_file} ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}
            COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}"
            VERBATIM
        )
        list(APPEND PYTHON_TEST_DEPS ${CMAKE_CURRENT_BINARY_DIR}/${py_file})
    endforeach ()

    add_custom_target(copy_python_tests ALL DEPENDS ${PYTHON_TEST_DEPS})
endif ()

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
