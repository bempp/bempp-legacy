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
lookup_python_package(cython REQUIRED PATH "${EXTERNAL_ROOT}/python")

cmake_test(pytest_build NOEXEC)
cmake_test(pytest_tests)

set_tests_properties(
    cmake_test_pytest_tests cmake_test_pytest_build
    PROPERTIES LABELS "python;pytest"
)
set_tests_properties(cmake_test_pytest_tests PROPERTIES
    DEPENDS cmake_test_pytest_build)
