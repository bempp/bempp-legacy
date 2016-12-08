include(PythonPackageLookup)
include(EnvironmentScript)

set(LOCAL_PYTHON_EXECUTABLE "${PROJECT_BINARY_DIR}/localpython.sh")
create_environment_script(
    EXECUTABLE "${PYTHON_EXECUTABLE}"
    PATH "${LOCAL_PYTHON_EXECUTABLE}"
    PYTHON
)
add_to_python_path("${EXTERNAL_ROOT}/python")

lookup_python_package(pytest REQUIRED PATH "${EXTERNAL_ROOT}/python")
find_program(PYTEST_EXEC py.test PATH "${EXTERNAL_ROOT}/python" REQUIRED)

cmake_test(pure_build NOEXEC --build-target install)
cmake_test(pure_tests)
set_tests_properties(cmake_test_pure_tests PROPERTIES
    DEPENDS cmake_test_pure_build)

cmake_test(extension_build NOEXEC --build-target install)
set(python_exec "${CMAKE_CURRENT_BINARY_DIR}/extension_tester.sh")
add_test(NAME python_extension
    COMMAND ${python_exec} ${PYTEST_EXEC} test_extension.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
set_tests_properties(python_extension PROPERTIES
    DEPENDS cmake_test_extension_build)
cmake_test(extension_files)
set_tests_properties(cmake_test_extension_files PROPERTIES
    DEPENDS cmake_test_extension_build)


cmake_test(cython_build NOEXEC --build-target install)
cmake_test(cython_files)
set_tests_properties(cmake_test_cython_files PROPERTIES
    DEPENDS cmake_test_cython_build)
set(python_exec "${CMAKE_CURRENT_BINARY_DIR}/cython_tester.sh")
add_test(NAME cython_extension
    COMMAND ${python_exec} ${PYTEST_EXEC} test_cython.py
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
set_tests_properties(cython_extension PROPERTIES
    DEPENDS cmake_test_cython_build)

set_tests_properties(
    cmake_test_pure_tests cmake_test_pure_build
    cmake_test_extension_build python_extension cmake_test_extension_files
    cmake_test_cython_build cmake_test_cython_files cython_extension
    PROPERTIES LABELS "python"
)
