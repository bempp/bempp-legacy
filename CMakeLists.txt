cmake_minimum_required(VERSION 2.8)
project(bempp)
enable_language(Fortran) # for BLAS

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -O3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DFIBER_CHECK_ARRAY_BOUNDS")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DFIBER_CHECK_ARRAY_BOUNDS")

# Module Path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Macros needed
include(BemppOptions)
include(BemppFindDependencies)

# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")


# Main library
add_subdirectory(lib)

# Python Support
add_subdirectory(python)

# Examples
add_subdirectory(examples) 

# Scratch
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scratch/CMakeLists.txt")
  add_subdirectory(scratch)
endif ()

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
