cmake_minimum_required(VERSION 2.6)
project(bempp)

# Options (can be modified by user)
option(WITH_TESTS "Compile unit tests (can be run with 'make test')" ON)

# Main library
add_subdirectory(lib) 

# Examples
add_subdirectory(examples) 

# Tests
if (WITH_TESTS)
  add_subdirectory(tests)
endif (WITH_TESTS)

# Uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
