# Figures out where to install python stuff

add_custom_target(cython-headers)
add_custom_target(cython)
add_subdirectory(bempp)
add_dependencies(cython cython-headers)

if(WITH_TESTS)
  file(WRITE ${PROJECT_BINARY_DIR}/bin/bempp_python_tests
    "${PROJECT_BINARY_DIR}/localpython.sh -c \"import bempp.api; bempp.api.test()\" ")
  execute_process(COMMAND chmod +x ${PROJECT_BINARY_DIR}/bin/bempp_python_tests)
  add_test(
    NAME python_tests
    COMMAND ${PROJECT_BINARY_DIR}/bin/bempp_python_tests
  )
endif()
