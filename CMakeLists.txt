cmake_minimum_required(VERSION 2.8)
project(bempp CXX)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")

# Tells cmake where to find bem++ specific files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Download and installs the cookoff in the build tree
# From here on we can use cookoff recipes
include(LookUp-GreatCMakeCookOff)

# Makes this a c++11 project
include(AddCPP11Flags)

# Macros needed
include(BemppOptions)
include(PackageLookUp)
include(BemppFindDependencies)

# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")


# Main library
add_subdirectory(lib)
#
# # Python Support
# add_subdirectory(python)
#
# Examples
add_subdirectory(examples)

# Scratch
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scratch/CMakeLists.txt")
  add_subdirectory(scratch)
endif ()

# Tests
if (WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif (WITH_TESTS)

# Uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
