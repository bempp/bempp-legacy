file(GLOB_RECURSE PYTHON_TEST_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach (py_file ${PYTHON_TEST_SOURCES})
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
    COMMAND ${CMAKE_COMMAND} copy -E ${CMAKE_CURRENT_SOURCE_DIR}/${py_file} ${CMAKE_CURRENT_BINARY_DIR}/${py_file}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}
    # TODO: remove the project's root directory from the printed path.
    COMMENT "Copying file ${CMAKE_CURRENT_SOURCE_DIR}/${py_file}"
    VERBATIM
)
endforeach ()

foreach(py_file ${PYTHON_TEST_SOURCES})
    set(PYTHON_TEST_DEPS "${PYTHON_TEST_DEPS};${CMAKE_CURRENT_BINARY_DIR}/${py_file}")
endforeach()
message(${PYTHON_TEST_DEPS})
add_custom_target(copy_python_tests ALL DEPENDS ${PYTHON_TEST_DEPS})

