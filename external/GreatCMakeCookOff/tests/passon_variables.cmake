find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
include(PassonVariables)

set(thisvar_zero "${CMAKE_CURRENT_BINARY_DIR}/hello" CACHE PATH "something")
set(thisvar_two True CACHE BOOL "something" FORCE)
set(othervar_one 2 CACHE STRING "something" FORCE)
set(othervar_two 2 CACHE INTERNAL "something" FORCE)
set(othervarone 3 CACHE PATH "something" FORCE)
set(alist 42;this;that CACHE STRING "something" FORCE)
passon_variables(thispackage
  FILENAME "${CMAKE_CURRENT_BINARY_DIR}/thispackage.cmake"
  PUBLIC
  PATTERNS ".*var_.*" alist
  ALSOADD
    "set(alsoadded \"hello world\")\n"
    "set(alsoadded2 42 this that)\n"
)

configure_file( "${cookoff_path}/tests/passon_variables_test.in.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/test.cmake"
    @ONLY
)

execute_process(
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/test.cmake"
    RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "passon_variables test failed -- ${result}")
endif()
