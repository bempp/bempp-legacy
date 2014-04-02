cmake_minimum_required(VERSION 2.8)
project(bempp)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")

# Module Path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Macros needed
include(BemppOptions)
include(PackageLookUp)
include(BemppFindDependencies)
#
# # set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")
#
# # Main library
# add_subdirectory(lib)
#
# # Python Support
# add_subdirectory(python)
#
# # Examples
# add_subdirectory(examples)
#
# # Scratch
# if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scratch/CMakeLists.txt")
#   add_subdirectory(scratch)
# endif ()
#
# # Tests
# if (WITH_TESTS)
#   add_subdirectory(tests)
# endif (WITH_TESTS)
#
# # Uninstall target
# configure_file(
#   "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
#   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
#   IMMEDIATE @ONLY)
#
# add_custom_target(uninstall
#   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
